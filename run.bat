@echo off
REM AI-Physics Drug Optimization Platform - Windows Startup Script
REM Usage: run.bat         (default port 7860)
REM        run.bat 8080    (custom port)

SET SCRIPT_DIR=%~dp0
SET PORT=%1
IF "%PORT%"=="" SET PORT=7860

echo [INIT] Checking Python...
python --version >nul 2>&1
IF ERRORLEVEL 1 goto no_python

echo [INIT] Checking dependencies...
python -c "import flask, flask_cors, Bio, numpy, scipy" >nul 2>&1
IF ERRORLEVEL 1 goto install_deps
goto check_db

:no_python
echo [ERROR] Python not found. Please install Python 3.9+ from https://www.python.org
pause
exit /b 1

:install_deps
echo [INIT] Installing required packages...
pip install flask flask-cors biopython numpy scipy -q
IF ERRORLEVEL 1 goto install_failed
goto check_db

:install_failed
echo [ERROR] Failed to install packages.
echo         Try manually: pip install flask flask-cors biopython numpy scipy
pause
exit /b 1

:check_db
IF NOT EXIST "%SCRIPT_DIR%db\skempi.db" goto no_db
goto start

:no_db
echo [ERROR] Database not found: %SCRIPT_DIR%db\skempi.db
echo         Rebuild with: python db\build_database.py
pause
exit /b 1

:start
echo [INIT] Database OK
echo [INIT] Starting server on http://localhost:%PORT%
echo [INFO] Open your browser to http://localhost:%PORT%
echo [INFO] Press Ctrl+C to stop.
echo.

cd /d "%SCRIPT_DIR%"
SET PORT=%PORT%
python api\server.py
pause
