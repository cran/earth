@rem test.earthmain.gcc.bat: test 64 bit standalone earth.c with main()
@rem
@rem This tests the earth C code.  It does this: builds test.earthmain.exe
@rem (under gcc), runs it, and compares results to test.earthmain.out.save.
@rem
@rem You will need to tweak this file for your directories.
@rem
@rem You need to make R.lib first -- see instructions in gnuwin32/README.packages.

@echo test.earthmain.gcc.bat
@set CYGWIN=nodosfilewarning
@rem Init environment for GCC compiler, if necessary
@call D:\bin\milbo\rpath.bat

cp "C:/Program Files/R/R-4.3.2/bin/x64/R.dll" .
                                @if %errorlevel% neq 0 goto err
cp "C:/Program Files/R/R-4.3.2/bin/x64/Rblas.dll" .
                                @if %errorlevel% neq 0 goto err
cp "C:/Program Files/R/R-4.3.2/bin/x64/Riconv.dll" .
                                @if %errorlevel% neq 0 goto err
cp "C:/Program Files/R/R-4.3.2/bin/x64/Rgraphapp.dll" .
                                @if %errorlevel% neq 0 goto err

@rem @rem you may have to create Rdll_x64.lib and Rblas_x64.lib beforehand
@cp "../../.#/Rdll_x64.lib" R.lib
                                 @if %errorlevel% neq 0 goto err
@cp "../../.#/Rblas_x64.lib" Rblas.lib
                                 @if %errorlevel% neq 0 goto err

@rem TODO -USE_BLAS=0 else crashes in daxpy_ call in FindKnot
@rem TODO -Wno-stringop-overflow else earth.c:3301:warning: 'memset' exceeds maximum object size

gcc -DMAIN -DSTANDALONE -DUSE_BLAS=0 -Wall --pedantic -Wextra -O3 -std=gnu99 -m64^
 -Wno-stringop-overflow^
 -I"/a/r/ra/include" -I../../inst/slowtests^
 ../../src/earth.c^
 R.lib Rblas.lib^
 -o earthmain.gcc.exe
                                @if %errorlevel% neq 0 goto err
earthmain.gcc.exe >test.earthmain.gcc.out
                                @rem no errorlevel test, diff will do check for discrepancies
                                @rem @if %errorlevel% neq 0 goto err
mks.diff test.earthmain.gcc.out test.earthmain.gcc.out.save
                                @if %errorlevel% neq 0 goto err

@rm -f R.dll Rblas.dll R.lib Rblas.lib Riconv.dll Rgraphapp.dll R.lib Rblas.lib
@rm -f earthmain.gcc.* test.earthmain.gcc.out *.o
@exit /B 0

:err
@exit /B %errorlevel%
