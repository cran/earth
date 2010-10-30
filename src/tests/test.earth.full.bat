@rem test.earth.full.bat: this does a regression test of earth
@rem Stephen Milborrow Apr 2007 Petaluma

@echo === test.earth.full ==============================================
@"\PROGRA~2\R\R-2.11.1\bin\R.exe" CMD BATCH --quiet --vanilla test.earth.full.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.full.Rout:
@echo.
@tail test.earth.full.Rout
@echo.
@exit /B 1
:good1
diff test.earth.full.Rout test.earth.full.Rout.save 
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@exit /B %errorlevel%
:good2
@rem test.earth.full.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diff -w Rplots.ps ..\..\.#\test.earth.full.save.ps 
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B %errorlevel%
:good3
@rm -f test.earth.full.Rout
@rm -f Rplots.ps
@exit /B  0
