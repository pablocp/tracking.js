(function() {
  /**
   * EPnp utility.
   * @static
   * @constructor
   */
  tracking.EPnP = {};

  tracking.EPnP.prototype.initPoints = function(objectPoints, imagePoints) {
    var instance = this,
        numberOfCorrespondences = instance.numberOfCorrespondences,
        pws = instance.pws,
        us = instance.us;

    for(var i = 0; i < numberOfCorrespondences; i++)
    {
      pws[3 * i    ] = objectPoints[i].x;
      pws[3 * i + 1] = objectPoints[i].y;
      pws[3 * i + 2] = objectPoints[i].z;

      us[2 * i    ] = imagePoints[i].x*fu + uc;
      us[2 * i + 1] = imagePoints[i].y*fv + vc;
    }
  };

  tracking.EPnP.prototype.initCameraParameters = function(cameraMatrix) {
    var instance = this;

    instance.uc = cameraMatrix[0][2];
    instance.vc = cameraMatrix[1][2];
    instance.fu = cameraMatrix[0][0];
    instance.fv = cameraMatrix[1][1];
  };

  tracking.EPnP.prototype.init = function(objectPoints, imagePoints, cameraMatrix) {
    var instance = this,
        numberOfCorrespondences = objectPoints.length;

    initCameraParameters(cameraMatrix);

    instance.numberOfCorrespondences = numberOfCorrespondences;
    instance.pws = new Float64Array(3*numberOfCorrespondences);
    instance.us = new Float64Array(2*numberOfCorrespondences);

    initPoints(objectPoints, imagePoints);

    instance.alphas = new Float64Array(4*numberOfCorrespondences);
    instance.pcs = new Float64Array(3*numberOfCorrespondences);

    instance.max_nr = 0;
  };

  // Decompose a m x n matrix using SVD
  tracking.EPnP.prototype.svd = function(A, m, n, W, U, V) {
    var matrix = [], i, j;

    for (i = 0; i < m; i++) {
      matrix.push([]);
      for (j = 0; j < n; j++) {
        matrix[i].push(A[i*n+j]);
      }
    }
    
    var output = numeric.svd(matrix),
        w = output.S,
        u = output.U,
        v = output.V;

    if (W) {
      for (i = 0; i < w.length; i++) {
        W[i] = w[i];
      }
    }

    if (U) {
      for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
          U[i*m + j] = u[i][j];
        }
      }
    }

    if (V) {
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
          V[i*n + j] = v[i][j];
        }
      }
    }
  };

  tracking.EPnP.prototype.invertSquare = function(src, n, dst) {
    var matrix = [], i, j;

    for (i = 0; i < n; i++) {
      matrix.push([]);
      for (j = 0; j < n; j++) {
        matrix[i].push(src[i*n+j]);
      }
    }

    matrix = numeric.inv(matrix);

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        dst[i*n + j] = matrix[i][j];
      }
    }
  }

  tracking.EPnP.prototype.transposeSquare = function(A, n) {
    var i, j, temp;

    for (i = 1; i < n; i++) {
      for (j = 0; j < i; j++) {
        temp = A[i*n+j];
        A[i*n+j] = A[j*n+i];
        A[j*n+i] = temp;
      }
    }
  };

  // Solves a linear system Ax = b
  tracking.EPnP.prototype.solveLinearSystem = function(A, m, n, b, dst) {
    var leftSide = [], 
        rightSide = [],
        i, j;

    for (i = 0; i < m; i++) {
      leftSide.push([]);
      for (j = 0; j < n; j++) {
        leftSide[i].push(A[i*n+j]);
      }
      rightSide.push(b[i]);
    }

    var output = numeric.solve(leftSide, rightSide);

    for (i = 0; i < n; i++) {
      dst[i] = output[i];
    }
  };

  tracking.EPnP.prototype.chooseControlPoints = function() {
    var instance = this,
        cws = new Float64Array(4*3),
        numberOfCorrespondences = instance.numberOfCorrespondences,
        position = 0,
        pws = instance.pws,
        i,
        j;

    instance.cws = cws;

    // Take C0 as the reference points centroid:
    cws[0] = cws[1] = cws[2] = 0;

    for(i = 0; i < numberOfCorrespondences; i++) {
      for(j = 0; j < 3; j++) {
        cws[position] += pws[position];
        ++position;
      }
    }

    for(j = 0; j < 3; j++)
      cws[j] /= numberOfCorrespondences;

    // Take C1, C2, and C3 from PCA on the reference points:
    var PW0 = new Float64Array(numberOfCorrespondences*3),
        PW0tPW0 = new Float64Array(3 * 3),
        DC = new Float64Array(3),
        UCt = new Float64Array(3 * 3);

    for(i = 0; i < numberOfCorrespondences; i++) {
      for(j = 0; j < 3; j++) {
        PW0[3 * i + j] = pws[3 * i + j] - cws[0][j];
      }
    }

    instance.mulTransposed(PW0, PW0tPW0, numberOfCorrespondences, 3, 1);

    instance.svd(PW0tPW0, 3, 3, DC, UCt, 0);
    instance.transposeSquare(UCt, 3);

    for(i = 1; i < 4; i++) {
      var k = sqrt(DC[i - 1] / numberOfCorrespondences);
      for(j = 0; j < 3; j++) {
        cws[i*3 + j] = cws[j] + k * UCt[3 * (i - 1) + j];
      }
    }
  };

  // Calculates the product of a m x n matrix and its transposition.
  tracking.EPnP.prototype.mulTransposed = function(src, m, n, dst, order) {
    var i, j, k;
    if(order) {
      // dst = srct x src
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
          dst[i*n + j] = 0;
          for (k = 0; k < m; k++) {
            dst[i*n + j] += src[k*n+i]*src[k*n+j];
          }
        }
      }
    }
    else {
      // dst = src x srct
      for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
          dst[i*n + j] = 0;
          for (k = 0; k < n; k++) {
            dst[i*n + j] += src[i*n+k]*src[j*n+k]; 
          }
        }
      }
    }
  };

  tracking.EPnP.prototype.computeBarycentricCoordinates = function() {
    var instance = this,
        alphas = instance.alphas,
        cc = new Float64Array(3*3),
        ccInv = new Float64Array(3*3),
        cws = instance.cws,
        i,
        j;

      for(i = 0; i < 3; i++) {
        for(j = 1; j < 4; j++) {
          cc[3 * i + j - 1] = cws[j][i] - cws[0][i];
        }
      }

      instance.invertSquare(cc, 3, ccInv);

      for(i = 0; i < instance.numberOfCorrespondences; i++) {
        var pi = 3 * i,;
            a = 4 * i;

        for(int j = 0; j < 3; j++) {
          alphas[a + 1 + j] =
            ccInv[3 * j    ] * (pws[pi + 0] - cws[0]) +
            ccInv[3 * j + 1] * (pws[pi + 1] - cws[1]) +
            ccInv[3 * j + 2] * (pws[pi + 2] - cws[2]);
        }
        alphas[a + 0] = 1.0f - alphas[a + 1] - alphas[a + 2] - alphas[a + 3];
      }
    };

    tracking.EPnP.prototype.fillM = function(M, row, as, offset, u, v) {
      var instance = this,
          fu = instance.fu,
          fv = instance.fv,
          uc = instance.uc,
          vc = instance.vc,
          m1 = row * 12,
          m2 = m1 + 12;

      for(int i = 0; i < 4; i++) {
        M[m1 + 3 * i    ] = as[offset + i] * fu;
        M[m1 + 3 * i + 1] = 0.0;
        M[m1 + 3 * i + 2] = as[offset + i] * (uc - u);

        M[m2 + 3 * i    ] = 0.0;
        M[m2 + 3 * i + 1] = as[offset + i] * fv;
        M[m2 + 3 * i + 2] = as[offset + i] * (vc - v);
      }
    };

    tracking.EPnP.prototype.computeCCS = function(betas, ut) {
      var instance = this,
          betas = instance.betas;
          ccs = new Float64Array(4 * 3),
          i,
          j,
          k;

      instance.ccs = ccs;

      for(i = 0; i < 4; i++) {
        ccs[i*3] = ccs[i*3 + 1] = ccs[i*3 + 2] = 0.0f;
      }

      for(i = 0; i < 4; i++) {
        var v = 12 * (11 - i);
        for(j = 0; j < 4; j++){
          for(k = 0; k < 3; k++){
            ccs[j*3 + k] += betas[i] * ut[v + 3 * j + k];
          }
        }
      }
    };

    tracking.EPnP.prototype.computePCS = function() {
      var instance = this,
          alphas = instance.alphas,
          ccs = instance.ccs,
          pcs = instance.pcs,
          i;

      for(i = 0; i < instance.numberOfCorrespondences; i++) {
        var a = 4 * i,
            pc = 3 * i;

        for(j = 0; j < 3; j++)
          pcs[pc + j] = alphas[a + 0] * ccs[0 * 3 * i + j] + 
                        alphas[a + 1] * ccs[1 * 3 * i + j] + 
                        alphas[a + 2] * ccs[2 * 3 * i + j] + 
                        alphas[a + 3] * ccs[3 * 3 * i + j];
      }
    };

    tracking.EPnP.prototype.computePose = function(R, t) {
      var instance = this,
          numberOfCorrespondences = instance.numberOfCorrespondences,
          i;

      instance.chooseControlPoints();
      instance.computeBarycentricCoordinates();

      var M = new Float64Array(2 * numberOfCorrespondences * 12);

      for(i = 0; i < numberOfCorrespondences; i++) {
        instance.fillM(M, 2 * i, alphas[0], 4 * i, us[2 * i], us[2 * i + 1]);
      }

      var MtM = new Float64Array(12*12),
          D = new Float64Array(12),
          Ut = new Float64Array(12*12);

      instance.mulTransposed(M, MtM, 2*numberOfCorrespondences, 12, 1);

      instance.svd(MtM, 12, 12, D, Ut, 0);
      instance.transposeSquare(Ut, 12);

      var L_6x10 = new Float64Array(6 * 10),
          Rho = new Float64Array(6);

      instance.computeL6x10(Ut, L_6x10);
      instance.computeRho(Rho);

      var Betas = [new Float64Array(4), new Float64Array(4), new Float64Array(4), new Float64Array(4)],
          repErrors = new Float64Array(4),
          Rs = [new Float64Array(3,3), new Float64Array(3,3), new Float64Array(3,3), new Float64Array(3,3)],
          ts = [new Float64Array(3), new Float64Array(3), new Float64Array(3), new Float64Array(3)];

      findBetasApprox1(L_6x10, Rho, Betas[1]);
      //gaussNewton(L_6x10, Rho, Betas[1]);
      repErrors[1] = computeRAndT(Ut, Betas[1], Rs[1], ts[1]);

      findBetasApprox2(L_6x10, Rho, Betas[2]);
      //gaussNewton(L_6x10, Rho, Betas[2]);
      repErrors[2] = computeRAndT(Ut, Betas[2], Rs[2], ts[2]);

      findBetasApprox3(L_6x10, Rho, Betas[3]);
      //gaussNewton(L_6x10, Rho, Betas[3]);
      repErrors[3] = computeRAndT(Ut, Betas[3], Rs[3], ts[3]);

      var N = 1;
      if (repErrors[2] < repErrors[1]) N = 2;
      if (repErrors[3] < repErrors[N]) N = 3;

      copyRAndT(Rs[N], Ts[n], R, T);
    };

    tracking.EPnP.prototype.copyRAndT = function(Rsrc, Tsrc, Rdst, Tdst) {
      var i, j;

      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
          Rdst[3*i + j] = Rsrc[3*i + j];
        }
        Tdst[i] = Tsrc[i];
      }
    };

    tracking.EPnP.prototype.dist2 = function(p1, p1offset, p2, p2offset) {
      return
        (p1[p1offset+0] - p2[p2offset+0]) * (p1[p1offset+0] - p2[p2offset+0]) +
        (p1[p1offset+1] - p2[p2offset+1]) * (p1[p1offset+1] - p2[p2offset+1]) +
        (p1[p1offset+2] - p2[p2offset+2]) * (p1[p1offset+2] - p2[p2offset+2]);
    };

    tracking.EPnP.prototype.dot = function(v1, v1offset, v2, v2offset) {
      return v1[v1offset+0] * v2[v2offset+0] + 
             v1[v1offset+1] * v2[v2offset+1] + 
             v1[v1offset+2] * v2[v2offset+2];
    };

    tracking.EPnP.prototype.estimateRAndT = function(R, t) {
      var instance = this,
          numberOfCorrespondences = instance.numberOfCorrespondences,
          pc0 = new Float64Array(3),
          pcs = instance.pcs,
          pw0 = new Float64Array(3),
          pws = instance.pws,
          i,
          j;

      pc0[0] = pc0[1] = pc0[2] = 0.0;
      pw0[0] = pw0[1] = pw0[2] = 0.0;

      for(i = 0; i < numberOfCorrespondences; i++) {
        var pc = 3 * i,
            pw = 3 * i;

        for(j = 0; j < 3; j++) {
          pc0[j] += pcs[pc + j];
          pw0[j] += pws[pw + j];
        }
      }
      for(j = 0; j < 3; j++) {
        pc0[j] /= numberOfCorrespondences;
        pw0[j] /= numberOfCorrespondences;
      }

      var ABt = new Float64Array(3 * 3), 
          ABt_D = new Float64Array(3), 
          ABt_U = new Float64Array(3 * 3), 
          ABt_V = new Float64Array(3 * 3);

      for (i = 0; i < 9; i++) {
        ABt = 0;
      }

      for(i = 0; i < numberOfCorrespondences; i++) {
        var pc = 3 * i,
            pw = 3 * i;

        for(int j = 0; j < 3; j++) {
          ABt[3 * j    ] += (pcs[pc + j] - pc0[j]) * (pws[pw + 0] - pw0[0]);
          ABt[3 * j + 1] += (pcs[pc + j] - pc0[j]) * (pws[pw + 1] - pw0[1]);
          ABt[3 * j + 2] += (pcs[pc + j] - pc0[j]) * (pws[pw + 2] - pw0[2]);
        }
      }

      instance.svd(ABt, 3, 3, ABt_D, UBt_U, ABt_V);

      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          R[i*3 + j] = instance.dot(ABt_U, 3 * i, ABt_V, 3 * j);
        }
      }

      const double det =
        R[0*3+ 0] * R[1*3+ 1] * R[2*3+ 2] + R[0*3+ 1] * R[1*3+ 2] * R[2*3+ 0] + R[0*3+ 2] * R[1*3+ 0] * R[2*3+ 1] -
        R[0*3+ 2] * R[1*3+ 1] * R[2*3+ 0] - R[0*3+ 1] * R[1*3+ 0] * R[2*3+ 2] - R[0*3+ 0] * R[1*3+ 2] * R[2*3+ 1];

      if (det < 0) {
        R[2*3+ 0] = -R[2*3+ 0];
        R[2*3+ 1] = -R[2*3+ 1];
        R[2*3+ 2] = -R[2*3+ 2];
      }

      t[0] = pc0[0] - dot(R, 0*3, pw0, 0);
      t[1] = pc0[1] - dot(R, 1*3, pw0, 0);
      t[2] = pc0[2] - dot(R, 2*3, pw0, 0);
    };

    tracking.EPnP.prototype.solveForSign = function() {
      var instance = this,
          pcs = instance.pcs,
          ccs = instance.ccs,
          i, 
          j;

      if (pcs[2] < 0.0) {
        for(i = 0; i < 4; i++) {
          for(j = 0; j < 3; j++) {
            ccs[i*3 + j] = -ccs[i*3 + j];
          }
        }

        for(i = 0; i < instance.numberOfCorrespondences; i++) {
          pcs[3 * i    ] = -pcs[3 * i];
          pcs[3 * i + 1] = -pcs[3 * i + 1];
          pcs[3 * i + 2] = -pcs[3 * i + 2];
        }
      }
    };

    tracking.EPnP.prototype.computeRAndT = function(ut, betas, R, t) {
      var instance = this;
      
      instance.computeCCS(betas, ut);
      instance.computePCS();

      instance.solveForSign();

      instance.estimateRAndT(R, t);

      return instance.reprojectionError(R, t);
    };

    tracking.EPnP.prototype.reprojectionError = function(R, t) {
      var instance = this,
          pws = instance.pws,
          dot = instance.dot,
          us = instance.us,
          uc = instance.uc,
          vc = instance.vc,
          fu = instance.fu,
          fv = instance.fv,
          sum2 = 0.0,
          i;

      for(i = 0; i < instance.numberOfCorrespondences; i++) {
        var pw = 3 * i,
            Xc = dot(R, 0*3, pws, pw) + t[0],
            Yc = dot(R, 1*3, pws, pw) + t[1],
            inv_Zc = 1.0 / (dot(R, 2*3, pws, pw) + t[2]),
            ue = uc + fu * Xc * inv_Zc,
            ve = vc + fv * Yc * inv_Zc,
            u = us[2 * i], v = us[2 * i + 1];

        sum2 += sqrt( (u - ue) * (u - ue) + (v - ve) * (v - ve) );
      }

      return sum2 / instance.numberOfCorrespondences;
    };

    // betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
    // betas_approx_1 = [B11 B12     B13         B14]
    tracking.EPnP.prototype.findBetasApprox1 = function(L_6x10, Rho, betas) {
      var L_6x4 = new Float64Array(6 * 4),
          B4 = new Float64Array(4),
          i;

      for(i = 0; i < 6; i++) {
        L_6x4[i*4 + 0] = L_6x10[i*10 + 0];
        L_6x4[i*4 + 1] = L_6x10[i*10 + 1];
        L_6x4[i*4 + 2] = L_6x10[i*10 + 3];
        L_6x4[i*4 + 3] = L_6x10[i*10 + 6];
      }

      instance.solve(L_6x4, 6, 4, Rho, B4);

      if (B4[0] < 0) {
        betas[0] = sqrt(-B4[0]);
        betas[1] = -B4[1] / betas[0];
        betas[2] = -B4[2] / betas[0];
        betas[3] = -B4[3] / betas[0];
      } else {
        betas[0] = sqrt(B4[0]);
        betas[1] = B4[1] / betas[0];
        betas[2] = B4[2] / betas[0];
        betas[3] = B4[3] / betas[0];
      }
    };

    // betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
    // betas_approx_2 = [B11 B12 B22                            ]
    tracking.EPnP.prototype.findBetasApprox2 = function(L_6x10, Rho, betas) {
      var L_6x3 = new Float64Array(6 * 3), 
          B3 = new Float64Array(3),
          i;

      for(i = 0; i < 6; i++) {
        L_6x3[i*3 + 0] = L_6x10[i*10 + 0];
        L_6x3[i*3 + 1] = L_6x10[i*10 + 1];
        L_6x3[i*3 + 2] = L_6x10[i*10 + 2];
      }

      instance.solve(L_6x3, 6, 3, Rho, B3);

      if (B3[0] < 0) {
        betas[0] = sqrt(-B3[0]);
        betas[1] = (B3[2] < 0) ? sqrt(-B3[2]) : 0.0;
      } else {
        betas[0] = sqrt(B3[0]);
        betas[1] = (B3[2] > 0) ? sqrt(B3[2]) : 0.0;
      }

      if (B3[1] < 0) betas[0] = -betas[0];

      betas[2] = 0.0;
      betas[3] = 0.0;
    };

    // betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
    // betas_approx_3 = [B11 B12 B22 B13 B23                    ]
    tracking.EPnP.prototype.findBetasApprox3 = function(fL_6x10, Rho, betas) {
      var L_6x5 = new Float64Array(6 * 5), 
          B5 = new Float64Array(5),
          i;

      for(i = 0; i < 6; i++) {
        L_6x5[i*5 + 0] = L_6x10[i*10 + 0;
        L_6x5[i*5 + 1] = L_6x10[i*10 + 1;
        L_6x5[i*5 + 2] = L_6x10[i*10 + 2;
        L_6x5[i*5 + 3] = L_6x10[i*10 + 3;
        L_6x5[i*5 + 4] = L_6x10[i*10 + 4;
      }

      instance.solve(L_6x5, 6, 5, Rho, B5);

      if (B5[0] < 0) {
        betas[0] = sqrt(-B5[0]);
        betas[1] = (B5[2] < 0) ? sqrt(-B5[2]) : 0.0;
      } else {
        betas[0] = sqrt(B5[0]);
        betas[1] = (B5[2] > 0) ? sqrt(B5[2]) : 0.0;
      }
      if (B5[1] < 0) betas[0] = -betas[0];
      betas[2] = B5[3] / betas[0];
      betas[3] = 0.0;
    };

    tracking.EPnP.prototype.computeL6x10 = function(ut, l_6x10) {
      var v = new UInt8Array(4),
          i;

      v[0] = 12 * 11;
      v[1] = 12 * 10;
      v[2] = 12 *  9;
      v[3] = 12 *  8;

      var dv = [new Float64Array(6*3), new Float64Array(6*3), new Float64Array(6*3), new Float64Array(6*3)];

      for(i = 0; i < 4; i++) {
        var a = 0, b = 1;
        for(int j = 0; j < 6; j++) {
          dv[i][j*3 + 0] = ut[v[i] + 3 * a    ] - ut[v[i] + 3 * b    ];
          dv[i][j*3 + 1] = ut[v[i] + 3 * a + 1] - ut[v[i] + 3 * b + 1];
          dv[i][j*3 + 2] = ut[v[i] + 3 * a + 2] - ut[v[i] + 3 * b + 2];

          b++;
          if (b > 3) {
            a++;
            b = a + 1;
          }
        }
      }

      for(i = 0; i < 6; i++) {
        l_6x10[10*i + 0] =        dot(dv[0], i*3, dv[0], i*3);
        l_6x10[10*i + 1] = 2.0f * dot(dv[0], i*3, dv[1], i*3);
        l_6x10[10*i + 2] =        dot(dv[1], i*3, dv[1], i*3);
        l_6x10[10*i + 3] = 2.0f * dot(dv[0], i*3, dv[2], i*3);
        l_6x10[10*i + 4] = 2.0f * dot(dv[1], i*3, dv[2], i*3);
        l_6x10[10*i + 5] =        dot(dv[2], i*3, dv[2], i*3);
        l_6x10[10*i + 6] = 2.0f * dot(dv[0], i*3, dv[3], i*3);
        l_6x10[10*i + 7] = 2.0f * dot(dv[1], i*3, dv[3], i*3);
        l_6x10[10*i + 8] = 2.0f * dot(dv[2], i*3, dv[3], i*3);
        l_6x10[10*i + 9] =        dot(dv[3], i*3, dv[3], i*3);
      }
    };

    tracking.EPnP.prototype.computeRho = function(rho) {
      var cws = this.cws;

      rho[0] = dist2(cws, 0*3, cws, 1*3);
      rho[1] = dist2(cws, 0*3, cws, 2*3);
      rho[2] = dist2(cws, 0*3, cws, 3*3);
      rho[3] = dist2(cws, 1*3, cws, 2*3);
      rho[4] = dist2(cws, 1*3, cws, 3*3);
      rho[5] = dist2(cws, 2*3, cws, 3*3);
    };

    tracking.EPnP.solve = function(objectPoints, imagePoints, cameraMatrix) {
      var instance = this;
    };

}());
