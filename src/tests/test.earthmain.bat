@rem test.earthmain.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Apr 2007 Petaluma

@echo === test.earthmain ===
@cp "C:\Program Files\r\R-2.5.0\bin\R.dll" .
@cp "C:\Program Files\r\R-2.5.0\bin\Rblas.dll" .
@md Debug

@cl -nologo -DSTANDALONE -DMAIN -TP -Zi -W3 -MLd -I"C:\Program Files\r\R-2.5.0\include" -I. -FpDebug\vc60.PCH -Fo"Debug/" -c ..\earth.c
@if %errorlevel% neq 0 goto error:
@link -nologo -debug:full -out:earthmain.exe Debug\earth.obj \a1\r\work\src\gnuwin32\Rdll.lib \a1\r\work\bin\Rblas.lib
@if %errorlevel% neq 0 goto error:
@earthmain.exe > Debug\test.earthmain.out
@if %errorlevel% neq 0 goto error:

@echo === Following diff should give no output ===
diff test.earthmain.out.save Debug\test.earthmain.out
@if %errorlevel% neq 0 goto error:

@rm -f R.dll Rblas.dll earthmain.exe *.map *.ilk *.pdb
@rm -rf Debug
@exit /B 0

:error
@echo error: nmake returned errorlevel %errorlevel%
@exit /B %errorlevel%
