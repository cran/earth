@rem test.earthmain.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Apr 2007 Petaluma

@echo test.earthmain.vc.bat
@set CYGWIN=nodosfilewarning

@mks.cp "C:\bin\R400devdll\i386\R.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Rblas.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Riconv.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Rgraphapp.dll" .
                                @if %errorlevel% neq 0 goto error
@rem you may have to create R.lib and Rblas.lib beforehand
mks.cp "C:\bin\R400devdll\i386\R.lib" .
                                @if %errorlevel% neq 0 goto error
mks.cp "C:\bin\R400devdll\i386\Rblas.lib" .
                                @if %errorlevel% neq 0 goto error

@md Debug

@rem Note: Use Microsoft VC14 (Visual Studio 2015) 32 bit.
@rem (Other versions haven't been tested and may cause spurious errors.)
@rem
@rem To set up the environment for the call to "cl" below, invoke:
@rem    C:\Program Files\Microsoft Visual Studio 14.0\Common7\Tools\vsvars32.bat
@rem
@rem We use -W4 below (insteadof -W3) for lint-like warnings

cl -nologo -DSTANDALONE -DMAIN -TP -Zi -W3 -MDd -I"%ProgramFiles%\R\R-4.0.2"\src\include -I. -FpDebug\vc60.PCH -Fo"Debug/" -c ..\..\src\earth.c
                                @if %errorlevel% neq 0 goto error
@rem Added path below to disambiguate from rtools which (may 2020 R version 4.0.0)
"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\BIN\link.exe" -nologo -debug -out:earthmain.exe Debug\earth.obj R.lib Rblas.lib
                                @if %errorlevel% neq 0 goto error
earthmain.exe > Debug\test.earthmain.out
                                @rem no errorlevel test, diff will do check for discrepancies
                                @rem @if %errorlevel% neq 0 goto error
mks.diff Debug\test.earthmain.out test.earthmain.out.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll R.lib Rblas.lib iconv.dll Riconv.dll Rgraphapp.dll earthmain.exe *.map *.ilk *.pdb
@rm -rf Debug
@exit /B 0

:error
@exit /B %errorlevel%
