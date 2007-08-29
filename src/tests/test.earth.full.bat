@rem test.earth.full.bat: this does a regression test of earth
@rem Stephen Milborrow Apr 2007 Petaluma

@echo === test.earth.full ===
"\Program Files\r\R-2.5.1\bin\Rcmd.exe" BATCH --quiet --vanilla test.earth.full.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.full.Rout
@exit /B %errorlevel%
:good1
@echo === Following two diffs should give no output ===
diff -w test.earth.full.Rout.save test.earth.full.Rout
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@exit /B %errorlevel%
@rem test.earth.full.save.ps is too big to be included in the release
@rem so it is stored elsewhere
:good2
diff -w ..\..\.#\test.earth.full.save.ps Rplots.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B %errorlevel%
:good3
@rm -f test.earth.full.Rout
@rm -f Rplots.ps
@exit /B  0
