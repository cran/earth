@rem earth/inst/slowtests/make.bat

time /T
@echo.
@call test.earthmain.vc.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.earthmain.gcc.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.earthmain.clang.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.earthc.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.glm.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.big.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.full.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.cv.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.plotd.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.varmod.bat
                        @if %errorlevel% NEQ 0 goto error
@echo.
@call test.weights.bat
                        @if %errorlevel% NEQ 0 goto error
@goto done
:error
@echo ==== ERROR ====
@exit /B %errorlevel%
:done
@rm -f ../../src/earth_res.rc ../Makedeps
@rm -f test.*.pdf *.dll *.lib *.pdb
time /T
@exit /B  0
