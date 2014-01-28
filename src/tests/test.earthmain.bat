@rem test.earthmain.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Apr 2007 Petaluma

@echo === test.earthmain ===============================================

@set CYGWIN=nodosfilewarning
@cp "C:/Program Files/R/R-3.0.2/bin/i386/R.dll" .
@cp "C:/Program Files/R/R-3.0.2/bin/i386/Rblas.dll" .
@cp "C:/a/r/ra/src/gnuwin32/unicode/iconv.dll" .
@cp "C:/Program Files/R/R-3.0.2/bin/i386/Riconv.dll" .
@cp "C:/Program Files/R/R-3.0.2/bin/i386/Rgraphapp.dll" .
@cp "C:/Program Files/R/R-3.0.2/bin/i386/Rzlib.dll" .
@rem you may have to create Rdll.lib and Rblas.lib beforehand
@cp "../../.#/Rdll.lib" .
@cp "../../.#/Rblas.lib" .
@md Debug

cl -nologo -DSTANDALONE -DMAIN -TP -Zi -W3 -I"C:/a/r/ra/include" -I. -FpDebug\vc60.PCH -Fo"Debug/" -c ..\earth.c
@if %errorlevel% neq 0 goto error:
link -nologo -debug -out:earthmain.exe Debug\earth.obj Rdll.lib Rblas.lib
@if %errorlevel% neq 0 goto error:
earthmain.exe > Debug\test.earthmain.out
@if %errorlevel% neq 0 goto error:

@rem we use -w on mks.diff so it treats \r\n the same as \n
mks.diff -w Debug\test.earthmain.out test.earthmain.out.save
@if %errorlevel% neq 0 goto error:

@rm -f R.dll Rblas.dll Rdll.lib Rblas.lib iconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll earthmain.exe *.map *.ilk *.pdb
@rm -rf Debug
@exit /B 0

:error
@echo error: errorlevel %errorlevel%
@exit /B %errorlevel%
