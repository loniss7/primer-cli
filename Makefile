PYTHON ?= python3
PKG_DIR := app/primer_cli
SRC_DIR := $(PKG_DIR)/primer_cli
TEST_DIR := $(PKG_DIR)/tests

.PHONY: format lint test typecheck check precommit-install precommit-run

format:
	$(PYTHON) -m ruff check --fix $(SRC_DIR) $(TEST_DIR)
	$(PYTHON) -m ruff format $(SRC_DIR) $(TEST_DIR)
	$(PYTHON) -m black --line-length 100 $(SRC_DIR) $(TEST_DIR)

lint:
	$(PYTHON) -m ruff check $(SRC_DIR) $(TEST_DIR)
	$(PYTHON) -m black --check --line-length 100 $(SRC_DIR) $(TEST_DIR)

typecheck:
	cd $(PKG_DIR) && $(PYTHON) -m mypy primer_cli tests

test:
	cd $(PKG_DIR) && $(PYTHON) -m pytest

check: lint typecheck test

precommit-install:
	$(PYTHON) -m pre_commit install

precommit-run:
	$(PYTHON) -m pre_commit run --all-files
