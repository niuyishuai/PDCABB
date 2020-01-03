# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 12:21:18 2018

@author: Yi-Shuai Niu 
email: niuyishuai@sjtu.edu.cn
Shanghai Jiao Tong University, 
Copyright (c) since 2017 Yi-Shuai Niu. All rights reserved.
"""

import ctypes as dll
import time

class DCABB_PARAL:
    """ Python interface for Paralle DCABB algorithm for solving mixed-integer optimization\n
        作者: Yi-Shuai Niu\n
        Shanghai Jiao Tong University, All rights reserved\n    
    """
    
    def __init__(self):
        '初始化'
        # 加载DCABB的dll
        self._dcabb=dll.CDLL('DCABB_PARAL_DLL.dll')
        # 定义常用函数结构
        # closedll
        self.closedll=dll.windll.kernel32.FreeLibrary
        self.closedll.argtypes=[dll.c_void_p]
        # DCABB::Create
        self.ReadModel=self._dcabb.ReadModel
        self.ReadModel.argtypes=[dll.c_char_p]
        # DCABB::Optimize
        self.Optimize=self._dcabb.Optimize
        # DCABB::GetSolutionStatus
        self.GetSolutionStatus=self._dcabb.GetSolutionStatus
        self.GetSolutionStatus.restype=dll.c_int
        # DCABB::GetObjVal
        self.GetObjVal=self._dcabb.GetObjVal
        self.GetObjVal.restype=dll.c_double              
        # DCABB::GetNbVars
        self.GetNbVars=self._dcabb.GetNbVars
        self.GetNbVars.restype=dll.c_int
        # DCABB::GetNbIntVars
        self.GetNbIntVars=self._dcabb.GetNbIntVars
        self.GetNbIntVars.restype=dll.c_int
        # DCABB::GetVarType
        self.GetVarType=self._dcabb.GetVarType
        self.GetVarType.argtypes=[dll.c_int]
        self.GetVarType.restype=dll.c_int
        # DCABB::GetX
        self.GetX=self._dcabb.GetX
        self.GetX.argtypes=[dll.c_int]
        self.GetX.restype=dll.c_double
        # DCABB::GetVarName
        self.GetVarName=self._dcabb.GetVarName
        self.GetVarName.argtypes=[dll.c_int]
        self.GetVarName.restype=dll.c_char_p
        # DCABB::GetCPUTime
        self.GetCPUTime=self._dcabb.GetCPUTime
        self.GetCPUTime.restype=dll.c_double
        # DCABB::GetNbDCA
        self.GetNbDCA=self._dcabb.GetNbDCA
        self.GetNbDCA.restype=dll.c_int
        # DCABB::GetNbNodes
        self.GetNbNodes=self._dcabb.GetNbNodes
        self.GetNbNodes.restype=dll.c_int
        # DCABB::GetNbIter
        self.GetNbIter=self._dcabb.GetNbIter
        self.GetNbIter.restype=dll.c_int
        # DCABB::SetTolGap
        self.SetTolGap=self._dcabb.SetTolGap
        self.SetTolGap.argtypes=[dll.c_double]
        self.SetTolGap.restype=dll.c_double
        # DCABB::SetTolGapRestartDCA
        self.SetTolGapRestartDCA=self._dcabb.SetTolGapRestartDCA
        self.SetTolGapRestartDCA.argtypes=[dll.c_double]
        self.SetTolGapRestartDCA.restype=dll.c_double
        # DCABB::SetParalMode
        self.SetParalMode=self._dcabb.SetParalMode
        self.SetParalMode.argtypes=[dll.c_bool]

         
    def __del__(self):
        time.sleep(1)
        self.closedll(self._dcabb._handle)
        
    def GetVars(self):
        self.xopt=[]
        nvars=self.GetNbVars()
        for i in range(nvars):
            self.xopt.append((self.GetVarName(i),self.GetX(i)))
        return self.xopt
