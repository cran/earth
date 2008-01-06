@rem test.all.bat

@call test.earthmain.gcc.bat
                        @if errorlevel 1 goto error
@call test.earthmain.bat
                        @if errorlevel 1 goto error
@call test.earthc.bat
                        @if errorlevel 1 goto error
@call test.earth.full.bat
                        @if errorlevel 1 goto error
@call test.plotmo.bat
                        @if errorlevel 1 goto error
@call test.earth.big.bat
                        @if errorlevel 1 goto error
@goto done
:error
@echo ==== ERROR ====
:done
