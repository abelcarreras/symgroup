      SUBROUTINE dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        IMPLICIT NONE
        ! Arguments
        INTEGER, INTENT(IN) :: n       ! Ordre of A matrix
        INTEGER, INTENT(IN) :: nrhs    ! number of vectors of independent term
        INTEGER, INTENT(IN) :: lda     ! dimension of A
        INTEGER, INTENT(IN) :: ldb     ! dimension of B
        INTEGER, INTENT(OUT) :: info   ! error code
        INTEGER, INTENT(OUT) :: ipiv(n)! permutation vector
        DOUBLE PRECISION, INTENT(INOUT) :: A(lda,n)  ! coefficients matrix
        DOUBLE PRECISION, INTENT(INOUT) :: B(ldb,nrhs)! independent term vector (solution)

        ! local variables
        INTEGER :: i, j, k, kp
        DOUBLE PRECISION :: temp, maxval
        DOUBLE PRECISION, ALLOCATABLE :: L(:,:)

        ! Inicialitzation
        info = 0

        ! check arguments
        IF (n < 0) THEN
          info = -1
          RETURN
        ELSE IF (nrhs < 0) THEN
          info = -2
          RETURN
        ELSE IF (lda < MAX(1,n)) THEN
          info = -4
          RETURN
        ELSE IF (ldb < MAX(1,n)) THEN
          info = -7
          RETURN
        END IF

        ! Allocate memory
        ALLOCATE(L(n,n))
        L = 0.0d0

        ! LU factorization
        DO k = 1, n-1
          ! find pivot
          maxval = ABS(A(k,k))
          kp = k
          DO i = k+1, n
            IF (ABS(A(i,k)) > maxval) THEN
              maxval = ABS(A(i,k))
              kp = i
            END IF
          END DO
          ipiv(k) = kp

          ! switch rows if needed
          IF (kp /= k) THEN
            DO j = 1, n
              temp = A(k,j)
              A(k,j) = A(kp,j)
              A(kp,j) = temp
            END DO
            DO j = 1, nrhs
              temp = B(k,j)
              B(k,j) = B(kp,j)
              B(kp,j) = temp
            END DO
          END IF

          ! delete Gaussian
          IF (A(k,k) /= 0.0d0) THEN
            DO i = k+1, n
              A(i,k) = A(i,k)/A(k,k)
              L(i,k) = A(i,k)
              DO j = k+1, n
                A(i,j) = A(i,j) - A(i,k)*A(k,j)
              END DO
            END DO
          ELSE
            info = k
            RETURN
          END IF
        END DO

        ! solve lower triangularsystem Ly = b
        DO j = 1, nrhs
          DO i = 2, n
            DO k = 1, i-1
              B(i,j) = B(i,j) - L(i,k)*B(k,j)
            END DO
          END DO
        END DO

        ! solve upper triangularsystem Ux = y
        DO j = 1, nrhs
          DO i = n, 1, -1
            temp = B(i,j)
            DO k = i+1, n
              temp = temp - A(i,k)*B(k,j)
            END DO
            B(i,j) = temp/A(i,i)
          END DO
        END DO

        DEALLOCATE(L)
      END SUBROUTINE dgesv