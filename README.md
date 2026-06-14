# Primer CLI

`primer-cli` — CLI-пайплайн для подбора праймеров к бактериальным генам с обязательной production-проверкой специфичности через локальную BLAST DB.

## Что делает пайплайн

Базовый контур:

1. `fetch` — загрузка CDS-последовательностей.
2. `align` — множественное выравнивание через MAFFT.
3. `conserved` — поиск консервативных окон.
4. `predict` — генерация и ранжирование пар праймеров.

Production-контур добавляет:

1. сборку или валидацию локальной BLAST DB;
2. FASTA QC после `fetch`;
3. BLAST-проверку уникальных праймеров;
4. pair-specificity policy для фильтрации кандидатных пар.

## Поддерживаемые режимы конфигурации

Поддерживаются два формата YAML:

- single-gene: один `project.target_gene`, явные `runtime.work_dir/output_dir/reports_dir/downloads_dir`;
- multi-gene: общий `runtime.root_dir` и список `genes`, где у каждого гена свои `fetch`, `specificity_db`, `design`, `blast_specificity`.

Multi-gene режим нужен, когда вы хотите прогнать несколько генов одним запуском и при этом оставить для каждого свои таксоны и пороги специфичности.

## Multi-gene конфиг

Актуальный пример: [config/examples/gene.production.yaml](/home/sassy/edia/GenePrimerSelect/config/examples/gene.production.yaml)

В нём показан batch-конфиг для `gene1` и `gene2`.

Структура:

- `project.name`, `project.version` — общие metadata batch-запуска.
- `runtime.root_dir` — общий корень, под которым автоматически создаются:
  `work/<gene>`, `out/<gene>`, `reports/<gene>`, `blastdb/<gene>_specificity_panel`.
- `runtime.shared_downloads_dir` — общий каталог ZIP-архивов NCBI Datasets для всех генов batch-запуска.
- `runtime.shared_unpack_dir` — общий каталог распакованных NCBI Datasets для всех генов batch-запуска.
- `tools` — пути к `datasets`, `mafft`, `blastn`, `makeblastdb`, `blastdbcmd`.
- `genes[]` — отдельная настройка на каждый ген.

Для каждого элемента `genes[]`:

- `gene` — имя гена;
- `fetch.query` или `query` — запрос для `fetch`;
- `specificity_db` — источники данных для BLAST DB именно этого гена;
- `design` — параметры design-этапа;
- `blast_specificity` — policy и BLAST thresholds именно этого гена.

`out_prefix` и `subjects_file` в multi-gene конфиге можно не указывать: они вычисляются автоматически из `runtime.root_dir`.
Если у нескольких генов есть общие таксоны, указывайте `runtime.shared_downloads_dir` и `runtime.shared_unpack_dir`, а `specificity_db.out_prefix` оставляйте отдельным для каждого гена.
Если `blast_specificity.require_predicted_on_target_amplicon: false`, `specificity_db.local_fasta` может быть пустым: тогда пайплайн работает как off-target-only проверка без локального target reference.

## NCBI Datasets assembly limits

`design.max_sequences` and `specificity_db.ncbi_datasets.assembly_limits` control different stages:

- `design.max_sequences` limits how many fetched gene/CDS sequences move into the design pipeline.
- `specificity_db.ncbi_datasets.assembly_limits` limits how many genome assemblies are selected and downloaded per taxon for BLAST DB construction.

The BLAST DB assembly limit is configured per taxon role and can be disabled with `null`:

```yaml
specificity_db:
  ncbi_datasets:
    assembly_level:
      - complete
    assembly_source: "RefSeq"
    annotated_only: true
    exclude_atypical: true
    exclude_multi_isolate: true
    exclude_mag: true
    assembly_limits:
      target: 30
      near_target: 15
      background: 5
```

When a role-specific limit is set, `primer-cli` first runs `datasets summary genome taxon`, normalizes and deduplicates assembly accessions, stores the selected accession list, downloads only those assemblies via `datasets download genome accession`, and writes a per-taxon manifest next to the cached ZIP.

The cache fingerprint includes the taxon, role, assembly filters, and the role-specific limit. Changing any of these values creates a new cache entry instead of reusing an incompatible archive.

## Как запускать

Запустить все гены из multi-gene конфига:

```bash
primer-cli production run-batch --config config/examples/gene.production.yaml
```

Запустить только один ген из multi-gene конфига:

```bash
primer-cli production run --config config/examples/gene.production.yaml --gene gene1
```

Собрать BLAST DB только для одного гена из multi-gene конфига:

```bash
primer-cli blastdb build --config config/examples/gene.production.yaml --gene gene1
```

Проверить готовую BLAST DB:

```bash
primer-cli blastdb validate --db runs/multi_gene_example/blastdb/gene1_specificity_panel
primer-cli blastdb info --db runs/multi_gene_example/blastdb/gene1_specificity_panel
```

Single-gene конфиги по-прежнему поддерживаются:

```bash
primer-cli production run --config path/to/gene.yaml
primer-cli blastdb build --config path/to/gene.yaml
```

## Что происходит при batch-запуске

Для каждого гена пайплайн выполняет один и тот же production-контур:

1. проверяет или пересобирает BLAST DB;
2. получает входные FASTA через `fetch` или `runtime.test_data_dir`;
3. запускает FASTA QC;
4. делает alignment;
5. находит консервативные окна;
6. строит пары праймеров;
7. запускает BLAST specificity;
8. пишет результаты в отдельные каталоги этого гена.

После batch-запуска создаётся общий файл:

- `reports/batch_summary.json`

В нём перечислены успешные и неуспешные гены, а также пути к их output/report/BLAST DB.

## Структура результатов для multi-gene

При `runtime.root_dir: ../../runs/multi_gene_example` структура будет такой:

```text
runs/multi_gene_example/
  blastdb/
    gene1_specificity_panel*
    gene2_specificity_panel*
  downloads/
    gene1/
    gene2/
  out/
    gene1/
    gene2/
  reports/
    batch_summary.json
    gene1/
    gene2/
  work/
    gene1/
    gene2/
```

Типичные production-артефакты для каждого гена:

- `out/<gene>/top_primers.csv`
- `out/<gene>/top_primers.json`
- `out/<gene>/top_primers.txt`
- `out/<gene>/reports/blast_summary.json`
- `out/<gene>/reports/pair_specificity.tsv`
- `reports/<gene>/<gene>_fetch_qc.json`
- `reports/<gene>/report_<gene>.json`
- `reports/<gene>/subjects.tsv`

## Offline и тестовый контур

Если задан `runtime.test_data_dir`, stage `fetch` не обращается к интернету, а ищет локальные FASTA:

- `<gene>_raw.fasta`
- `<gene>.raw.fasta`
- `<gene>.fasta`
- `<gene>.fa`
- `raw.fasta`

Это позволяет прогонять production-контур на синтетических данных.

Пример offline-конфига: [config/examples/test_gene.synthetic_offline.production.yaml](/home/sassy/edia/GenePrimerSelect/config/examples/test_gene.synthetic_offline.production.yaml)

## Логирование

- любая команда принимает `--log-level` и `--log-file`;
- `production run` и `blastdb build` пишут логи в `reports/<gene>/logs/` или `reports_dir/logs/`;
- `production run-batch` пишет общий batch-log в `reports/logs/`, а также отдельные логи каждого гена.

Пример:

```bash
primer-cli production run-batch --config config/examples/gene.production.yaml --log-level DEBUG
```

## Команды

- `fetch`
- `align`
- `conserved`
- `predict`
- `run`
- `blastdb build`
- `blastdb validate`
- `blastdb info`
- `production run`
- `production run-batch`

## Требования

- Python 3.10+
- MAFFT
- для production-режима: `datasets`, `makeblastdb`, `blastdbcmd`, `blastn`
- доступ в интернет, если не используется `runtime.test_data_dir`

## Документация

- [docs/blast_specificity.md](/home/sassy/edia/GenePrimerSelect/docs/blast_specificity.md)
