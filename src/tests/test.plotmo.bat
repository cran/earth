@rem test.plotmo.bat: this does a regression test of plotmo
@rem Stephen Milborrow Apr 2007 Petaluma

@echo === test.plotmo ==================================================
@"\PROGRA~1\R\R-2.12.0\bin\R.exe" CMD BATCH --quiet --vanilla test.plotmo.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.plotmo.Rout:
@echo.
@tail test.plotmo.Rout
@echo.
@exit /B 1
:good1
diff test.plotmo.Rout test.plotmo.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test.plotmo.save.ps
@exit /B 1
:good2
@rem test.plotmo.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test.plotmo.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.plotmo.Rout
@rm -f Rplots.ps
@exit /B 0
