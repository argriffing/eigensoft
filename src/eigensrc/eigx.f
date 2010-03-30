                SUBROUTINE EIGX(MAT, W, N)
                INTEGER N 
                DOUBLE PRECISION MAT(*), W(*)  
                CHARACTER  JOBZ, UPLO
                INTEGER    INFO, LDZ, I
                DOUBLE PRECISION Z(N,N), WORK(3*N)
                JOBZ = 'N'
                UPLO = 'L' 
                LDZ = N
 13             FORMAT(F9.3)
                CALL DSPEV  (JOBZ, UPLO, N, MAT, W, Z, LDZ, WORK, INFO )
                IF (INFO.EQ.0) RETURN
                WRITE(6,11) INFO
 11             FORMAT('INFO:',I6)  
                END

                SUBROUTINE EIGXV(MAT, W, Z, N)
                INTEGER N 
                DOUBLE PRECISION MAT(*), W(*)  , Z(N, N)
                CHARACTER  JOBZ, UPLO
                INTEGER    INFO, LDZ, I
                DOUBLE PRECISION WORK(3*N)
                JOBZ = 'V'
                UPLO = 'L' 
                LDZ = N
 13             FORMAT(F9.3)
                CALL DSPEV  (JOBZ, UPLO, N, MAT, W, Z, LDZ, WORK, INFO )
                IF (INFO.EQ.0) RETURN
                WRITE(6,11) INFO
 11             FORMAT('INFO:',I6)  
                END

                SUBROUTINE    HELLO
                 WRITE(6,101)
 101             FORMAT('hello world from fortran')  
                END
