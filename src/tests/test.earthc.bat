@rem test.earthc.bat:
@rem
@rem This tests the earth C code.  It does this: builds test.earthc.exe
@rem (under MS VC 6.0), runs it, and compares results to test.earthc.out.save
@rem You need to make Rdll.lib first -- see instructions in gnuwin32/README.packages
@rem You will need to tweak this file and test.earthc.mak for your directories
@rem
@rem Stephen Milborrow Mar 2007 Forden, Wales

@echo === test.earthc ==================================================
@cp "C:/Program Files/R/R-2.15.0/bin/i386/R.dll" .
@cp "C:/Program Files/R/R-2.15.0/bin/i386/Rblas.dll" .
@cp "C:/a/r/ra/src/gnuwin32/unicode/iconv.dll" .
@cp "C:/Program Files/R/R-2.15.0/bin/i386/Riconv.dll" .
@cp "C:/Program Files/R/R-2.15.0/bin/i386/Rgraphapp.dll" .
@cp "C:/Program Files/R/R-2.15.0/bin/i386/Rzlib.dll" .
@rem you may have to create Rdll.lib and Rblas.lib beforehand
@cp "../../.#/Rdll.lib" .
@cp "../../.#/Rblas.lib" .
@md Debug
@md Release

@nmake -nologo CFG=Release -f test.earthc.mak %1 %2 %3

@rem The Debug build gives slightly different output in lower decimal places (TODO why?)
@rem The advantage of using Debug is that memory leaks are reported.
@rem It is much slower though.
@rem nmake -nologo CFG=Debug -f test.earthc.mak %1 %2 %3

@if %errorlevel% equ 0 goto good:
@echo error: errorlevel %errorlevel%
@exit /B %errorlevel%
:good
@rm -f R.dll Rblas.dll Rdll.lib Rblas.lib iconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll test.earthc.main.exe test.earthc.main.map test.earthc.main.ilk *.pdb
@rm -rf Debug
@rm -rf Release
