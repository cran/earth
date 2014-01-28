@rem test.earth.full.bat: this does a regression test of earth
@rem Stephen Milborrow Apr 2007 Petaluma

@echo === test.earth.full ==============================================
@"\PROGRA~1\R\R-3.0.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.earth.full.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.full.Rout:
@echo.
@tail test.earth.full.Rout
@echo.
@exit /B 1
:good1
mks.diff test.earth.full.Rout test.earth.full.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.earth.full.save.ps
@exit /B 1
:good2
@rem test.earth.full.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.earth.full.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.earth.full.Rout
@rm -f Rplots.ps
@exit /B  0
