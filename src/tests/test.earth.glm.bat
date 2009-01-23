@rem test.earth.glm.bat

@echo === test.earth.glm ===============================================
@"\a\r\ra\bin\Rcmd.exe" BATCH --quiet --vanilla test.earth.glm.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.glm.Rout
@exit /B %errorlevel%
:good1
diff test.earth.glm.Rout test.earth.glm.Rout.save 
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@exit /B %errorlevel%
:good2
@rem test.earth.glm.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diff -w Rplots.ps ..\..\.#\test.earth.glm.save.ps 
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B %errorlevel%
:good3
@rm -f test.earth.glm.Rout
@rm -f Rplots.ps
@exit /B  0
