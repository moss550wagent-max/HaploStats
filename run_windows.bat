@echo off
echo Checking and installing HaploStats dependencies...
python -m pip install --upgrade pip
python -m pip install fastapi uvicorn pandas openpyxl pydantic requests numpy

echo.
echo Starting HaploStats Clinical Engine on http://localhost:8000 ...
python -m uvicorn scripts.api:app --host 127.0.0.1 --port 8000 --reload
pause
