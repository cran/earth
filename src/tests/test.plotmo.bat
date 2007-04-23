@rem test.plotmo.bat: this does a regression test of plotmo
@rem Stephen Milborrow Apr 2007 Petaluma

"c:\Program Files\r\R-2.4.1\bin\Rcmd.exe" BATCH --quiet --vanilla test.plotmo.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.plotmo.Rout
@exit /B %errorlevel%
:good1
@echo The following two diffs should give no output
diff test.plotmo.Rout.save test.plotmo.Rout
@if %errorlevel% equ 0 goto good2:
@echo error: files are different
@exit /B %errorlevel%
@rem test.plotmo.save.ps is too big to be included in the release
@rem so it is stored elsewhere
:good2
diff ..\..\.#\test.plotmo.save.ps Rplots.ps 
@if %errorlevel% equ 0 goto good3:
@echo error: files are different
@exit /B %errorlevel%
:good3
@rm -f test.plotmo.Rout
@rm -f Rplots.ps
@exit /B 0
