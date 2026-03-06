# Contributing to microbiome-disease-pipeline

Thank you for your interest in contributing!

## How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Development Setup

```bash
git clone https://github.com/jmatos10/microbiome-disease-pipeline.git
cd microbiome-disease-pipeline
pip install -r requirements.txt
pip install -e ".[dev]"
```

## Code Style

- Follow PEP 8 (max line length: 120)
- Use type hints for all function signatures
- Add docstrings (Google style) to all public functions
- Write tests for new features

## Running Tests

```bash
pytest tests/ -v
flake8 src/ tests/ --max-line-length=120
```

## Reporting Bugs

Open an issue with steps to reproduce, expected vs actual behavior, and environment details.
