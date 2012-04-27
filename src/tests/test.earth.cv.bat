@rem test.earth.cv.bat: tests earth cross validation
@rem Stephen Milborrow Nov 2008 Gardens

@echo === test.earth.cv ==============================================
@"\PROGRA~1\R\R-2.15.0\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.earth.cv.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.cv.Rout:
@echo.
@tail test.earth.cv.Rout
@echo.
@exit /B 1
:good1
diff test.earth.cv.Rout test.earth.cv.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.earth.cv.save.ps
@exit /B 1
@rem test.earth.cv.save.ps is too big to be included in the release
@rem so it is stored elsewhere
:good2
diffps Rplots.ps ..\..\.#\test-reference\test.earth.cv.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.earth.cv.Rout
@rm -f Rplots.ps
@exit /B  0
