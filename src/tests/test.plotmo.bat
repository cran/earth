@rem test.plotmo.bat: this does a regression test of plotmo
@rem Stephen Milborrow Apr 2007 Petaluma

@echo === test.plotmo ==================================================
@"\a\r\ra\bin\Rcmd.exe" BATCH --quiet --vanilla test.plotmo.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.plotmo.Rout
@exit /B %errorlevel%
:good1
diff test.plotmo.Rout test.plotmo.Rout.save 
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@exit /B %errorlevel%
:good2
@rem test.plotmo.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diff -w Rplots.ps ..\..\.#\test.plotmo.save.ps 
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B %errorlevel%
:good3
@rm -f test.plotmo.Rout
@rm -f Rplots.ps
@exit /B 0
