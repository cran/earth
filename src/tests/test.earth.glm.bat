@rem test.earth.glm.bat

@echo === test.earth.glm ===============================================
@"\PROGRA~1\R\R-2.15.0alpha\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.earth.glm.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.glm.Rout:
@echo.
@tail test.earth.glm.Rout
@echo.
@exit /B 1
:good1
diff test.earth.glm.Rout test.earth.glm.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.earth.glm.save.ps
@exit /B 1
:good2
@rem test.earth.glm.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.earth.glm.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.earth.glm.Rout
@rm -f Rplots.ps
@exit /B  0
