# Theoretical Calculations

Набор Python-модулей для расчётов, связанных с рентгеновской лабораторной практикой: свойства материалов, оптические константы, коэффициенты рассеяния и вспомогательные базы данных.

## Что есть в проекте

- `src/materials` - работа со свойствами веществ и элементов. Модуль читает локальные SQLite-базы, умеет получать свойства соединений и при необходимости подтягивать данные из Materials Project.
- `src/coefficients` - расчёт оптических констант `delta`, `beta`, `n`, `e1`, `e2`, `mu`, а также производных параметров для зеркальных пар absorber/spacer.
- `src/coefficients/plotting` - функции построения графиков Plotly и Streamlit-приложение для интерактивной работы.
- `src/coefficients/calc` - расчёты optical jump, максимального отражения, resolving power и penetration depth.
- `data` - локальные данные проекта: SQLite-базы compounds/elements/PICS/atomic factors и таблицы `f1f2` для элементов.
- `src/definitions.py` - централизованные настройки путей и кредов через `pydantic-settings` и `.env`.

В репозитории также есть старые/исследовательские директории вроде `photoionization`, `atomic_profiles`, `global_modules`, `thickness_estimation`, `wagner_plot`. Они сохранены как рабочие материалы и legacy-инструменты, но основная актуальная пакетная структура находится в `src`.

## Streamlit-приложение

Приложение запускается из корня проекта:

```bash
streamlit run streamlit_app.py
```

В приложении доступны:

- зависимость оптических констант от энергии или длины волны;
- график зависимости `δ` от `β` при выбранной энергии;
- расчёт параметров зеркала для одной или нескольких пар материалов.

## Настройки

Секреты и переопределения путей читаются из `.env` в корне проекта. Пример:

```env
MAPI_KEY=
F1F2_DIR=
```

`MAPI_KEY` нужен только для обращения к Materials Project, если локальная база не содержит нужного соединения.

## Установка

Основной вариант для локальной разработки:

```bash
poetry install
poetry run streamlit run streamlit_app.py
```

Для Streamlit Cloud добавлены:

- `streamlit_app.py` в корне проекта;
- `requirements.txt`;
- `runtime.txt`;
- `.streamlit/config.toml`.
