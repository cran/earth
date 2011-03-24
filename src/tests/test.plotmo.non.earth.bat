@rem test.plotmo.non.earth.bat: test plotmo on non-earth models
@rem Stephen Milborrow, Basley KwaZulu-Natal Mar 2011

@echo === test.plotmo.non.earth ==================================================
@"\PROGRA~1\R\R-2.12.0\bin\R.exe" CMD BATCH --quiet --vanilla test.plotmo.non.earth.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.plotmo.non.earth.Rout:
@echo.
@tail test.plotmo.non.earth.Rout
@echo.
@exit /B 1
:good1
diff test.plotmo.non.earth.Rout test.plotmo.non.earth.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test.plotmo.non.earth.save.ps
@exit /B 1
:good2
@rem test.plotmo.non.earth.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test.plotmo.non.earth.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.plotmo.non.earth.Rout
@rm -f Rplots.ps
@exit /B 0
