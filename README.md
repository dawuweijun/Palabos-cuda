Palabos With Cuda Support
=============
*   [Overview](#overview)

*   [一，Palabos的基本架构：](#Palabos的基本架构)
    *   [1,Palabos的数据结构;](#Some Data Structures Of Palabos)
    *   [2,Palabos的并行方案;](#)
        *   [a.OpenMPI和OpenMP;](#)
        *   [b.CoProcessor;](#)
        *   [c.Why we need Internal GPU-based Palabos ?](#)
    *   [3.Palabos架构的优缺点;](#)
    
*   [二，GPU并行编程有哪些？](#)
    * [1, OpenACC](#)
    * [2, OpenCL](#)
    * [3, CUDA](#)
    * [4, Why we choose CUDA?](#)
    
*   [RoadMap of CUDA support for Palabos](#)

*   [Current Status:](#)
    *   [探讨技术路线;](#)
    *   [Work On this Document](#)
    *   [Waiting for Cuda 6;](#)
    
    Cuda 6 即将提供对Unified Memory(统一内存地址) 的支持，这将极大降低异构编程的复杂性，
为了降低编码的复杂性，决定对Cuda支持的工作要等到2014年Cuda发布之后进行。
