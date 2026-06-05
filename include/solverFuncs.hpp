#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include <sstream>
#include "HYPRE_struct_ls.h"
#include "HYPRE_struct_mv.h"
#include "HYPRE_krylov.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include <mpi.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>

using namespace std;

//*****************************************************************
//**                                                             **
//**           void solver_SMG()                                 **
//**                                                             **
//*****************************************************************
void solver_SMG(std::vector<std::vector<double> > aP
    ,std::vector<std::vector<double> > aE
    ,std::vector<std::vector<double> > aW
    ,std::vector<std::vector<double> > aN
    ,std::vector<std::vector<double> > aS
    ,std::vector<std::vector<double> > b
    ,std::vector<std::vector<int> > i_bl
    ,std::vector<std::vector<int> > j_bl
    ,int n_bl,int maxITR,double error_cnv,int icon_msg
    ,std::vector<std::vector<double> > &solution
){

    HYPRE_Init();
    // 格子の定義
    HYPRE_StructGrid grid;
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);

    for (int iblock=0;iblock<n_bl;iblock++){
        HYPRE_Int ilower[2] = {i_bl[iblock][0], j_bl[iblock][0]};
        HYPRE_Int iupper[2] = {i_bl[iblock][1], j_bl[iblock][1]};
        HYPRE_StructGridSetExtents(grid, ilower, iupper);
    }

    HYPRE_StructGridAssemble(grid);

    // ステンシルの定義
    HYPRE_StructStencil stencil;
    HYPRE_StructStencilCreate(2, 5, &stencil);

    HYPRE_Int offsets[5][2] = {{ 0,  0},
                         { 1,  0},
                         {-1,  0},
                         { 0,  1},
                         { 0, -1}};
    for (int i = 0; i < 5; i++) {
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }

    // 行列の作成
    HYPRE_StructMatrix A;
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    // 格子点ごとに異なる値を設定
    HYPRE_Int stencil_indices[5] = {0,1,2,3,4};
    //std::cout << ilower1[0] << ", "<<iupper1[0] << ","<<ilower1[1] << ", "<<iupper1[1] << std::endl;
    for (int iblock=0;iblock<n_bl;iblock++){
        for (int i=i_bl[iblock][0];i<=i_bl[iblock][1];i++){
            for (int j=j_bl[iblock][0];j<=j_bl[iblock][1];j++){
                HYPRE_Int index[2] = {i, j};
                double values[5];
                values[0] = aP[i][j];  // 中央
                values[1] = -aE[i][j];  // 右
                values[2] = -aW[i][j];  // 左
                values[3] = -aN[i][j];  // 上
                values[4] = -aS[i][j];  // 下
                HYPRE_StructMatrixSetValues(A, index, 5, stencil_indices, values);
            }
        }
    }

    HYPRE_StructMatrixAssemble(A);

    // ベクトルの作成
    HYPRE_StructVector bvec;
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &bvec);
    HYPRE_StructVectorInitialize(bvec);

    for (int iblock=0;iblock<n_bl;iblock++){
        for (int i=i_bl[iblock][0];i<=i_bl[iblock][1];i++){
            for (int j=j_bl[iblock][0];j<=j_bl[iblock][1];j++){
                HYPRE_Int index[2] = {i, j};
                double rhs = b[i][j];  // 右辺を1に設定
                HYPRE_StructVectorSetValues(bvec, index, rhs);
            }
        }
    }

    HYPRE_StructVectorAssemble(bvec);

    HYPRE_StructVector xvec;
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &xvec);
    HYPRE_StructVectorInitialize(xvec);
    HYPRE_StructVectorAssemble(xvec);

    // SMGソルバの作成
    HYPRE_StructSolver solver;
    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);

    // SMGソルバの設定
    HYPRE_StructSMGSetTol(solver, error_cnv);
    HYPRE_StructSMGSetMaxIter(solver, maxITR);
    if(icon_msg==1){
        HYPRE_StructSMGSetPrintLevel(solver, 2);
        HYPRE_StructSMGSetLogging(solver, 1);
    }

    // ソルバのセットアップと解法
    HYPRE_StructSMGSetup(solver, A, bvec, xvec);
    clock_t time1 = clock();
    HYPRE_StructSMGSolve(solver, A, bvec, xvec);
    clock_t time2 = clock();
    double time_real = (double)(time2-time1)/CLOCKS_PER_SEC;
    
    // イタレーション数の出力
    if(icon_msg==1){
        std::cout << "solver time = " << time_real << " sec" << std::endl;
        HYPRE_Int num_iterations;
        HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
        std::cout << "Number of Iterations: " << num_iterations << std::endl;
    
        // 結果の出力
        double final_res_norm;
        HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        std::cout << "Final Relative Residual Norm: " << final_res_norm << std::endl;
    }

    for (int iblock=0;iblock<n_bl;iblock++){
        for (int i=i_bl[iblock][0];i<=i_bl[iblock][1];i++){
            for (int j=j_bl[iblock][0];j<=j_bl[iblock][1];j++){
                HYPRE_Int index[2] = {i, j};
                double value;
                HYPRE_StructVectorGetValues(xvec, index, &value);
                solution[i][j] = value;
            }
        }
    }

    // クリーンアップ
    HYPRE_StructSMGDestroy(solver);
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(bvec);
    HYPRE_StructVectorDestroy(xvec);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructGridDestroy(grid);
    HYPRE_Finalize();
}

//*****************************************************************
//**                                                             **
//**           void solver_SOR()                                 **
//**                                                             **
//*****************************************************************
void solver_SOR(std::vector<std::vector<double> > aP
    ,std::vector<std::vector<double> > aE
    ,std::vector<std::vector<double> > aW
    ,std::vector<std::vector<double> > aN
    ,std::vector<std::vector<double> > aS
    ,std::vector<std::vector<double> > b
    ,std::vector<std::vector<int> > i_bl
    ,std::vector<std::vector<int> > j_bl
    ,int n_bl, double alpha_SOR,int maxITR,double error_cnv,int icon_msg
    ,std::vector<std::vector<double> > &solution
){

    int ncount = 0;
    double error = 0.0;

    //clock_t time1 = clock();
    do{
        ncount = ncount+1;
        error = 0.0;

        for (int iblock=0;iblock<n_bl;iblock++){
            for (int i=i_bl[iblock][0];i<=i_bl[iblock][1];i++){
                for (int j=j_bl[iblock][0];j<=j_bl[iblock][1];j++){
                    //std::cout << " i = "<< i << " j = " << j << std::endl;
                    double solution_tmp = (aE[i][j]*solution[i+1][j] + aW[i][j]*solution[i-1][j] 
                        + aN[i][j]*solution[i][j+1] + aS[i][j]*solution[i][j-1] + b[i][j])/(aP[i][j]+1e-100);
                    double solution_new = (1.0-alpha_SOR)*solution[i][j] +  alpha_SOR*solution_tmp;
                    double error_tmp = fabs((solution_new-solution[i][j])/(solution_new+solution[i][j]+1e-6))*2.0;
                    if(error_tmp > error){
                        error = error_tmp;
                    }
                    solution[i][j] = solution_new;
                }
            }
        }

        //if(icon_msg==1){
        //    std::cout <<"ncount = " << ncount<< " error_max = " << error << std::endl;
        //}
        if(ncount > maxITR) break;

    }while(error > error_cnv);

    if(icon_msg==1){
        std::cout <<"ncount = " << ncount<< " error_max = " << error << std::endl;
    }

}

//*****************************************************************
//**                                                             **
//**           void solver_LU_9p()                               **
//**                                                             **
//*****************************************************************
void solver_LU_9p(std::vector<std::vector<double> > aP
    ,std::vector<std::vector<double> > aE
    ,std::vector<std::vector<double> > aW
    ,std::vector<std::vector<double> > aN
    ,std::vector<std::vector<double> > aS
    ,std::vector<std::vector<double> > aEN
    ,std::vector<std::vector<double> > aES
    ,std::vector<std::vector<double> > aWN
    ,std::vector<std::vector<double> > aWS
    ,std::vector<std::vector<double> > b
    ,std::vector<std::vector<int> > i_bl
    ,std::vector<std::vector<int> > j_bl
    ,int n_bl, int icon_msg
    ,std::vector<std::vector<double> > &solution
){

    int ni = solution.size()-2;
    int nj = solution[0].size()-2;

    //対応番号表
    std::vector<int> ik;
    std::vector<int> jk;
    std::vector<std::vector<int> > kpos(ni+2,std::vector<int>(nj+2,-1));

    /*
    for (int i=1;i<=ni;i++){
        for (int j=1;j<=nj;j++){
            int k = (i-1)*nj + (j-1);
            kpos[i][j] = k;
            ik.push_back(i);
            jk.push_back(j);
            //std::cout << "i = " << i << " j = " << j<< " k = " << k << std::endl;
        }
    }
    */

    int k=0;
    for (int iblock=0;iblock<n_bl;iblock++){
        for (int i=i_bl[iblock][0];i<=i_bl[iblock][1];i++){
            for (int j=j_bl[iblock][0];j<=j_bl[iblock][1];j++){
                kpos[i][j] = k;
                ik.push_back(i);
                jk.push_back(j);
                k++;
            }
        }
    }
    
    int nk = ik.size();
    //std::cout << "nk = " << nk << std::endl;
    //abort();

    Eigen::SparseMatrix<double> A(nk, nk);
    A.reserve(Eigen::VectorXi::Constant(nk,20)); //ここの数字を変えて帯域幅を確保する
    Eigen::VectorXd bv(nk);
    Eigen::VectorXd xv(nk);
        
    for (int k=0;k<nk;k++){
        int i = ik[k];
        int j = jk[k];
        //std::cout << "i = "<< i<< " j = "<< j << " k = "<< k  << std::endl;

        double kE = kpos[i+1][j];
        double kW = kpos[i-1][j];
        double kN = kpos[i][j+1];
        double kS = kpos[i][j-1];

        double kEN = kpos[i+1][j+1];
        double kES = kpos[i+1][j-1];
        double kWN = kpos[i-1][j+1];
        double kWS = kpos[i-1][j-1];
        //std::cout << " nk =" << nk << " kE = "<< kE<< " kW = "<< kW 
        //    << " kN = "<< kN << " kS = "<< kS  << std::endl;
  
        A.insert(k, k)       =  aP[i][j];
        bv[k]                =   b[i][j];
        if(kE!=-1){
            A.insert(k, kE)  = -aE[i][j];
        }
        if(kW!=-1){
            A.insert(k, kW)  = -aW[i][j];
        }
        if(kN!=-1){
            A.insert(k, kN)  = -aN[i][j];
        }
        if(kS!=-1){
            A.insert(k, kS)  = -aS[i][j];
        }
        if(kEN!=-1){
            A.insert(k, kEN) = -aEN[i][j];
        }
        if(kES!=-1){
            A.insert(k, kES) = -aES[i][j];
        }
        if(kWN!=-1){
            A.insert(k, kWN) = -aWN[i][j];
        }
        if(kWS!=-1){
            A.insert(k, kWS) = -aWS[i][j];
        }
        
        /*
        A.insert(k, k)       =  aP[i][j];
        bv[k]                +=   b[i][j];
        if(kE!=-1){
            A.insert(k, kE)  =  - aE[i][j];
        }
        
        if(kW!=-1){
            A.coeffRef(k, kW)  -= aW[i][j];
        }
        if(kN!=-1){
            A.coeffRef(k, kN)  -= aN[i][j];
        }
        if(kS!=-1){
            A.coeffRef(k, kS)  -= aS[i][j];
        }
        if(kEN!=-1){
            A.coeffRef(k, kEN) -= aEN[i][j];
        }
        if(kES!=-1){
            A.coeffRef(k, kES) -= aES[i][j];
        }
        if(kWN!=-1){
            A.coeffRef(k, kWN) -= aWN[i][j];
        }
        if(kWS!=-1){
            A.coeffRef(k, kWS) -= aWS[i][j];
        }
        */

    }


    // LU分解でAx = bを解く
    Eigen::SparseLU<Eigen::SparseMatrix<double > > solver;
    //std::cout << "solver start!" << std::endl;
    solver.compute(A);
    xv = solver.solve(bv);

    if(icon_msg == 1){
        // 実際の誤差を計算
        Eigen::VectorXd residual = A*xv- bv;
        //std::cout << residual<< std::endl;
        double actual_error = residual.norm()/(bv.norm()+1e-100);
        std::cout << "Actual error: " << actual_error << std::endl;
    }

    //Ex-結果を戻す
    for (int k=0;k<nk;k++){
        int i = ik[k];
        int j = jk[k];
        solution[i][j] = xv[k];
    }

}
