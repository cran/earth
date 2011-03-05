@rem test.all.bat

time /T
@echo.
@call test.earthmain.gcc.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@echo.
@call test.earthmain.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.earthc.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.earth.glm.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.plotmo.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.earth.big.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.earth.full.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.earth.cv.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.plotd.bat
                        @if %errorlevel% NEQ 0 goto error
@goto done
:error
@echo ==== ERROR ====
:done
time /T
@rm -f ../earth_res.rc ../Makedeps
@rm -f test.*.pdf
