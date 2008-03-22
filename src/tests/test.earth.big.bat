@rem test.earth.big.bat: This tests earth on a biggish model
@rem This is the test mentioned in the earth man page "Big Models" section
@rem Stephen Milborrow Mar 2008 Durban

@echo === test.earth.big ===
"\a\r\ra\bin\Rcmd.exe" BATCH --quiet --vanilla test.earth.big.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.big.Rout
@exit /B %errorlevel%
:good1
@echo === Following two diffs should give no output ===
diff test.earth.big.Rout.save test.earth.big.Rout
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@exit /B %errorlevel%
@rem test.earth.big.save.ps is too big to be included in the release
@rem so it is stored elsewhere
:good2
diff -w ..\..\.#\test.earth.big.save.ps Rplots.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B %errorlevel%
:good3
@rm -f test.earth.big.Rout
@rm -f Rplots.ps
@exit /B  0
