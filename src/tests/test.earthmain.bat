@rem test.earthmain.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Apr 2007 Petaluma

@echo === test.earthmain ===============================================

@set CYGWIN=nodosfilewarning
@cp "C:/a/r/ra/bin/R.dll" .
@cp "C:/a/r/ra/bin/Rblas.dll" .
@cp "C:/a/r/ra/bin/iconv.dll" .
@cp "C:/a/r/ra/bin/Riconv.dll" .
@cp "C:/a/r/ra/bin/Rgraphapp.dll" .
@cp "C:/a/r/ra/bin/Rzlib.dll" .
@md Debug

cl -nologo -DSTANDALONE -DMAIN -TP -Zi -W3 -MLd -I"C:/a/r/ra/include" -I. -FpDebug\vc60.PCH -Fo"Debug/" -c ..\earth.c
@if %errorlevel% neq 0 goto error:
link -nologo -debug:full -out:earthmain.exe Debug\earth.obj "C:\a\r\ra\bin\Rdll.lib" "C:\a\r\ra\bin\Rblas.lib"
@if %errorlevel% neq 0 goto error:
earthmain.exe > Debug\test.earthmain.out
@if %errorlevel% neq 0 goto error:

@rem we use -w on diff so it treats \r\n the same as \n
diff -w Debug\test.earthmain.out test.earthmain.out.save
@if %errorlevel% neq 0 goto error:

@rm -f R.dll Rblas.dll iconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll earthmain.exe *.map *.ilk *.pdb
@rm -rf Debug
@exit /B 0

:error
@echo error: errorlevel %errorlevel%
@exit /B %errorlevel%
