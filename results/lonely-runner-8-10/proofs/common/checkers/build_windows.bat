@echo off
rem Builds checkers\bin\gate_lrat_check.exe with the Visual Studio 2022 clang-cl.
set "PATH=%PATH%;C:\Program Files (x86)\Microsoft Visual Studio\Installer"
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\Llvm\x64\bin\clang-cl.exe" /nologo /O2 /EHsc /std:c++17 /D_CRT_SECURE_NO_WARNINGS gate_lrat_check.cpp /Fe:bin\gate_lrat_check.exe /Fo:bin\gate_lrat_check.obj
echo build exit %ERRORLEVEL%
