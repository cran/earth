@rem test.all.bat

@echo.
@call test.earthmain.gcc.bat
                        @if errorlevel 1 goto error
@echo.
@call test.earthmain.bat
                        @if errorlevel 1 goto error
@echo.
@call test.earthc.bat
                        @if errorlevel 1 goto error
@echo.
@call test.earth.glm.bat
                        @if errorlevel 1 goto error
@echo.
@call test.earth.full.bat
                        @if errorlevel 1 goto error
@echo.
@call test.plotmo.bat
                        @if errorlevel 1 goto error
@echo.
@call test.earth.big.bat
                        @if errorlevel 1 goto error

@goto done
:error
@echo ==== ERROR ====
:done

@rm -f test.*.pdf
