# Primer CLI

`primer-cli` — CLI-пайплайн для подбора праймеров по генам устойчивости и другим целевым генам на основе последовательностей из NCBI.

Пайплайн объединяет шаги:
1. `fetch` — загрузка CDS из NCBI E-utilities.
2. `align` — множественное выравнивание MAFFT.
3. `conserved` — поиск консервативных окон.
4. `predict` — генерация, фильтрация и ранжирование пар праймеров.

Команда `run` выполняет все шаги end-to-end.

## Возможности

- Полный запуск в одну команду: `primer-cli run`.
- Работа с одним или несколькими генами (`--genes "vanA,vanB,mcr-1"`).
- Экспорт результатов в `CSV`, `JSON`, текстовый отчёт.
- Гибкая параметризация термодинамики, покрытий и критериев 3'-конца.
- Отдельные команды для пошаговой отладки (`fetch`, `align`, `conserved`, `predict`).

## Структура репозитория

- `app/primer_cli/` — Python-пакет и исходный код CLI.
- `app/primer_cli/primer_cli/` — модули приложения.
- `app/primer_cli/tests/` — unit и integration тесты.
- `environment.yml` — окружение Conda (включая внешние bio-tools).
- `requirements.txt` — минимальная установка через pip (`-e ./app/primer_cli`).

## Требования

- Python `3.10+` (рекомендуется `3.11`).
- MAFFT в `PATH` (обязательно для выравнивания).
- Доступ в интернет для запросов к NCBI.
- Опционально: `blast`/`blast+` для расширенных проверок специфичности.

## Установка

### Вариант 1: Conda (рекомендуется)

```bash
conda env create -f environment.yml
conda activate biopy
```

### Вариант 2: venv + pip

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## Быстрый старт

```bash
source .venv/bin/activate
export NCBI_EMAIL="you@example.com"

primer-cli run \
  --genes "vanA" \
  --work-dir ./work \
  --output-dir ./out \
  --max-sequences 100
```

Минимально обязательные аргументы для `run`:
- `--genes`
- `--work-dir`
- `--output-dir`

Если `work-dir` и `output-dir` не существуют, они создаются автоматически.

## Команда run (полный пайплайн)

`primer-cli run` — это основной режим для базового пользователя.  
Одна команда выполняет весь цикл:
1. поиск и загрузка CDS из NCBI (`fetch`);
2. множественное выравнивание (`align`);
3. поиск консервативных окон (`conserved`);
4. подбор и ранжирование пар праймеров (`predict`).

Базовый шаблон:

```bash
primer-cli run \
  --genes "vanA" \
  --work-dir ./work \
  --output-dir ./out
```

### Ключевые параметры run

| Параметр | Обязательный | Что делает | Пример |
|---|---|---|---|
| `--genes` | Да | Ген или список генов через запятую | `"vanA"` / `"vanA,vanB,mcr-1"` |
| `--work-dir` | Да | Папка для промежуточных файлов (`raw.fasta`, `aligned.fasta`) | `./work` |
| `--output-dir` | Да | Папка с финальными результатами (`regions.json`, `top_primers.*`) | `./out` |
| `--max-sequences` | Нет | Сколько последовательностей скачать на ген (быстрее отладка на меньшем числе) | `50`, `100`, `300` |
| `--window-size` | Нет | Размер окна для поиска консервативных участков | `25` |
| `--top-quantile` | Нет | Порог отбора лучших консервативных окон (0..1] | `0.8` |
| `--top-n` | Нет | Сколько лучших пар праймеров сохранить в итогах | `20` |
| `--email` | Нет* | Email для NCBI (можно вместо этого `NCBI_EMAIL`) | `you@example.com` |
| `--mafft` | Нет | Путь к MAFFT (если не в `PATH`) | `/usr/bin/mafft` |
| `--validate-blast` | Нет | Включает дополнительную BLAST-проверку off-target | флаг |
| `--blast-db` | Нет** | Путь/имя BLAST нуклеотидной БД | `./data/blast_test_db/test_subjects` |

\* Практически рекомендуется всегда задавать `--email` или `NCBI_EMAIL`.
\** Обязателен только если включен `--validate-blast`.

### Быстрые примеры run

Один ген, стандартный запуск:

```bash
primer-cli run --genes "vanA" --work-dir ./work --output-dir ./out
```

Быстрый тестовый прогон (меньше данных):

```bash
primer-cli run --genes "vanA" --work-dir ./work --output-dir ./out --max-sequences 20
```

Несколько генов за один запуск:

```bash
primer-cli run --genes "vanA,vanB,mcr-1" --work-dir ./work --output-dir ./out --max-sequences 100
```

Запуск с включенной BLAST-валидацией:

```bash
primer-cli run \
  --genes "vanA" \
  --work-dir ./work \
  --output-dir ./out \
  --validate-blast \
  --blast-db ./data/blast_test_db/test_subjects
```

## Запуск по нескольким генам

```bash
primer-cli run \
  --genes "vanA,vanB,mcr-1" \
  --work-dir ./work \
  --output-dir ./out
```

Для multi-gene запуска результаты раскладываются по подпапкам:

```text
work/
  vanA/
    raw.fasta
    aligned.fasta
  vanB/
    raw.fasta
    aligned.fasta

out/
  vanA/
    regions.json
    top_primers.csv
    top_primers.json
    top_primers.txt
  vanB/
    ...
```

## Пошаговый режим

```bash
# 1) Скачивание последовательностей
primer-cli fetch --gene vanA --output ./work/raw.fasta --max-sequences 100

# 2) Выравнивание
primer-cli align --input ./work/raw.fasta --output ./work/aligned.fasta --mafft mafft

# 3) Консервативные окна
primer-cli conserved \
  --input ./work/aligned.fasta \
  --output ./out/regions.json \
  --window-size 25 \
  --top-quantile 0.8

# 4) Предсказание и ранжирование праймеров
primer-cli predict \
  --raw-fasta ./work/raw.fasta \
  --aligned-fasta ./work/aligned.fasta \
  --conserved-regions ./out/regions.json \
  --output-dir ./out
```

## Формат выходных файлов

В `output-dir` сохраняются:
- `regions.json` — найденные консервативные области.
- `top_primers.csv` — табличный рейтинг пар праймеров.
- `top_primers.json` — детализированный JSON с метриками.
- `top_primers.txt` — человекочитаемый отчёт.

Основные поля для интерпретации:
- `pair_coverage` — доля последовательностей, где пара проходит критерии.
- `tm_forward`, `tm_reverse` — температуры плавления (чем ближе друг к другу, тем лучше).
- `gc_forward`, `gc_reverse` — GC-состав праймеров.
- `amplicon_length` — длина ампликона.
- `final_score` — итоговый скоринг (0..100).
- `offtarget_summary` — статус проверки off-target.

## Логи и отладка

Запуск с подробными логами:

```bash
primer-cli --log-level DEBUG run --genes vanA --work-dir ./work --output-dir ./out
```

Если кажется, что процесс "завис":
- `fetch` может долго работать на больших выборках из NCBI;
- `align` (MAFFT) часто самый долгий этап;
- проверьте активность процессов `primer-cli` / `mafft`.

## Тесты

Из директории `app/primer_cli`:

```bash
pytest
```

## Частые проблемы

- `mafft: command not found`
  - Установите MAFFT и убедитесь, что бинарник доступен в `PATH`.

- Ошибки доступа к NCBI
  - Проверьте интернет, повторите запуск позже, задайте `--email` или `NCBI_EMAIL`.

- Путь указывает на файл вместо папки
  - Для `--work-dir` и `--output-dir` требуется директория, а не файл.

## Roadmap

- Расширенная автоматическая off-target проверка (BLAST) в основном пайплайне.
- Параллелизация отдельных стадий для ускорения multi-gene запуска.
- Расширение отчётности и визуализации качества праймеров.
