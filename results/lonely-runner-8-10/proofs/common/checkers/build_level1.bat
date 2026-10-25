@echo off
set "PATH=%PATH%;C:\Program Files (x86)\Microsoft Visual Studio\Installer"
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\Llvm\x64\bin\clang-cl.exe" /nologo /O2 /EHsc /std:c++17 /D_CRT_SECURE_NO_WARNINGS level1.cpp /Fe:bin\level1.exe /Fo:bin\level1.obj
echo build exit %ERRORLEVEL%
