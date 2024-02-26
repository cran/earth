@rem test.earthc.gcc.bat:
@rem
@rem This tests the earth C code.  It does this: builds test.earthc.exe
@rem (under gcc), runs it, and compares results to test.earthc.out.save.
@rem
@rem You will need to tweak this file for your directories.
@rem
@rem You need to make R.lib first -- see instructions in gnuwin32/README.packages.

@echo test.earthc.gcc.bat
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

gcc -DSTANDALONE -DUSE_BLAS=0 -Wall --pedantic -Wextra -O3 -std=gnu99 -m64^
 -Wno-stringop-overflow -Wno-unused-parameter^
 -I"/a/r/ra/include" -I../../inst/slowtests^
 ../../src/earth.c test.earthc.c^
 R.lib Rblas.lib^
 -o test.earthc.exe
                                @if %errorlevel% neq 0 goto err
test.earthc.exe >test.earthc.gcc.out
                                @if %errorlevel% neq 0 goto err
@rem we use -w on diff so it treats \r\n the same as \n
diff -w test.earthc.gcc.out test.earthc.gcc.out.save
                                @if %errorlevel% neq 0 goto err
@if %errorlevel% equ 0 goto good
:err
@echo error: errorlevel %errorlevel%
@exit /B %errorlevel%
:good
@rm -f R.dll Rblas.dll R.lib Rblas.lib Riconv.dll Rgraphapp.dll R.lib Rblas.lib
@rm -f test.earthc.exe test.earthc.gcc.out
@exit /B 0
