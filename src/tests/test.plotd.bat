@rem test.plotd.bat: This tests earth on a biggish model
@rem This is the test mentioned in the earth man page "Big Models" section
@rem Stephen Milborrow Mar 2008 Durban

@echo === test.plotd ===============================================
@"\a\r\ra\bin\Rcmd.exe" BATCH --quiet --vanilla test.plotd.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.plotd.Rout
@exit /B %errorlevel%
:good1
diff test.plotd.Rout test.plotd.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@exit /B %errorlevel%
:good2
@rem test.plotd.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diff -w Rplots.ps ..\..\.#\test.plotd.save.ps 
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B %errorlevel%
:good3
@rm -f test.plotd.Rout
@rm -f Rplots.ps
@exit /B  0