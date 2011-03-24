@rem test.plotmo.fac.bat: test factor plotting in plotmo. This also tests swapxy, xflip, and yflip
@rem Stephen Milborrow, Berea Mar 2011

@echo === test.plotmo.fac ==================================================
@"\PROGRA~1\R\R-2.12.0\bin\R.exe" CMD BATCH --quiet --vanilla test.plotmo.fac.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.plotmo.fac.Rout:
@echo.
@tail test.plotmo.fac.Rout
@echo.
@exit /B 1
:good1
diff test.plotmo.fac.Rout test.plotmo.fac.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test.plotmo.fac.save.ps
@exit /B 1
:good2
@rem test.plotmo.fac.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test.plotmo.fac.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.plotmo.fac.Rout
@rm -f Rplots.ps
@exit /B 0
