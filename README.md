# python-3DOF-robot-inverse-Dynamics-symbolic-operation
The process of solving the inverse dynamics of the puma type robot through symbolic operation
# Robot Dynamics Toolkit

이 프로젝트는 **로봇 동역학/기구학 계산**을 위한 Python 기반 툴킷입니다.  
CasADi와 Sympy를 활용하여 **관성 행렬(Inertia Tensor) 계산**,  
**Denavit-Hartenberg(DH) 변환 행렬 구성**, **야코비안 기반 질량행렬(Mass Matrix) 산출**을 수행할 수 있습니다.  

---

## 📦 설치 방법

### 1. Python 환경 준비
- Python 3.8 이상 권장

### 2. 필수 라이브러리 설치
```bash
pip install casadi
