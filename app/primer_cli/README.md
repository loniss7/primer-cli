# primer-cli

CLI-инструмент для подбора праймеров из MSA по пайплайну:

`NCBI fetch -> MAFFT align -> conserved windows -> primer prediction`

## Установка

```bash
pip install -e .
```

Или из корня репозитория:

```bash
pip install -e ./app/primer_cli
```

## Команды

```bash
primer-cli --help
primer-cli fetch --help
primer-cli align --help
primer-cli conserved --help
primer-cli predict --help
primer-cli run --help
```

## Полный запуск

```bash
export NCBI_EMAIL="you@example.com"

primer-cli run \
  --genes "vanA" \
  --work-dir ./work \
  --output-dir ./out \
  --max-sequences 100
```

Для нескольких генов:

```bash
primer-cli run --genes "vanA,vanB,mcr-1" --work-dir ./work --output-dir ./out
```

Результаты сохраняются в `output-dir` (`top_primers.csv`, `top_primers.json`, `top_primers.txt`) и промежуточные данные в `work-dir`.

## Разработка

```bash
pytest
```

## Примечания

- Требуется Python `3.10+`.
- Для стадии выравнивания нужен `mafft` в `PATH`.
- Для загрузки из NCBI нужен доступ в интернет.
