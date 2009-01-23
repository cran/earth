@rem test.earth.cv.bat: tests earth cross validation
@rem Stephen Milborrow Nov 2008 Gardens

@echo === test.earth.cv ==============================================
@"\a\r\ra\bin\Rcmd.exe" BATCH --quiet --vanilla test.earth.cv.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.cv.Rout
@exit /B %errorlevel%
:good1
diff test.earth.cv.Rout test.earth.cv.Rout.save 
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@exit /B %errorlevel%
@rem test.earth.cv.save.ps is too big to be included in the release
@rem so it is stored elsewhere
:good2
diff -w Rplots.ps ..\..\.#\test.earth.cv.save.ps 
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B %errorlevel%
:good3
@rm -f test.earth.cv.Rout
@rm -f Rplots.ps
@exit /B  0
