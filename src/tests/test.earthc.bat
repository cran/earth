@rem test.earthc.bat:
@rem 
@rem This tests the earth C code.  It does this: builds test.earthc.main.exe 
@rem (under MS VC 6.0), runs it, and compares results to test.earthc.ref
@rem You need to make Rdll.lib first -- see instructions in gnuwin32\README.packages
@rem You will need to tweak this file and test.earthc.mak for your directories
@rem 
@rem Stephen Milborrow Mar 2007 Forden, Wales

@set R_HOME \a1\r\work
@cp /a1/r/work/bin/R.dll .
@cp /a1/r/work/bin/Rblas.dll .
@md Debug
@md Release
nmake -nologo CFG=Release -f test.earthc.mak %1 %2 %3
@rem debug build gives slightly different output 
@rem nmake -nologo CFG=Debug -f test.earthc.mak %1 %2 %3
@if %errorlevel% equ 0 goto good:
@echo error: nmake returned errorlevel %errorlevel%
@exit /B %errorlevel%
:good
@rm -f R.dll Rblas.dll test.earthc.main.exe test.earthc.main.map test.earthc.main.ilk *.pdb
@rm -rf Debug
@rm -rf Release
