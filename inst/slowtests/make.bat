@rem earth/inst/slowtests/make.bat

@call test.earthmain.gcc.bat
                        @if %errorlevel% NEQ 0 goto error
@rem TODO Feb 2020: clang gives link errors if used with ms c compiler on the path
@rem @call test.earthmain.clang.bat
                   @if %errorlevel% NEQ 0 goto error
@call test.earthmain.vc.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.earthc.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.mods.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.incorrect.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.big.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.weights.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.expand.bpairs.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.bpairs.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.full.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.glm.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.allowedfunc.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.cv.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.pmethod.cv.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.varmod.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.varmod.mgcv.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.plotd.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.offset.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.ordinal.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.multresp.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.emma.bat
                        @if %errorlevel% NEQ 0 goto error
@rem TODO Sep 2020: commented out because varying results per run with R 4.0.2
@rem @call test.mem.bat
@rem                    @if %errorlevel% NEQ 0 goto error
@goto done
:error
@echo ==== ERROR ====
@exit /B %errorlevel%
:done
@rm -f ../../src/earth_res.rc ../Makedeps
@rm -f test.*.pdf *.dll *.lib *.pdb
@exit /B  0
