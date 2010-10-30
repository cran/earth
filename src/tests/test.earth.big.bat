@rem test.earth.big.bat: This tests earth on a biggish model
@rem This is the test mentioned in the earth man page "Big Models" section
@rem Stephen Milborrow Mar 2008 Durban

@echo === test.earth.big ===============================================
@"\PROGRA~2\R\R-2.11.1\bin\R.exe" CMD BATCH --quiet --vanilla test.earth.big.R
@if %errorlevel% equ 0 goto good1:
@echo error: R returned errorlevel %errorlevel%, see test.earth.big.Rout:
@echo.
@tail test.earth.big.Rout
@echo.
@exit /B 1
:good1
diff test.earth.big.Rout test.earth.big.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@exit /B %errorlevel%
:good2
@rem test.earth.big.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diff -w Rplots.ps ..\..\.#\test.earth.big.save.ps 
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B %errorlevel%
:good3
@rm -f test.earth.big.Rout
@rm -f Rplots.ps
@exit /B  0
