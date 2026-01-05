@echo off
echo ========================================
echo 银杏内酯B 分子对接分析
echo ========================================
echo.

REM 检查Docker是否运行
docker version >nul 2>&1
if %errorlevel% neq 0 (
    echo [错误] Docker未运行，请先启动Docker Desktop
    pause
    exit /b 1
)

REM 检查镜像是否存在
docker images | findstr autodock-vina >nul 2>&1
if %errorlevel% neq 0 (
    echo [提示] 首次运行，正在构建Docker镜像...
    docker build -t autodock-vina .
)

echo [1/3] 检查文件...
if not exist "ligand" mkdir ligand
if not exist "receptor" mkdir receptor
if not exist "results" mkdir results
if not exist "config" mkdir config

echo [2/3] 运行分子对接...
docker run --rm -v "%CD%:/docking" autodock-vina python /docking/molecular_docking.py

echo.
echo [3/3] 完成！
echo 结果保存在: results\docking_results.csv
echo.
pause
