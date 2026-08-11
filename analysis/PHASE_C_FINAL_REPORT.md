# Phase C 完了報告

## C3.3c 監査結果

### Audit 1: τ→0 limit ✅ PASS
- Γ_F^exact (cut2) = 4.000000
- C_dyn^(τ) → 4.000000 as τ→0
- factor of 2 のバグを修正（gamma_F.py の gamma_exact 関数）

### Audit 2: Frame covariance ✅ PASS
- Lab frame と F frame で C_F(0) が一致（Δ = 2.22e-16）

### Audit 3: D_F definition consistency ✅ PASS
- C3.3b の D_F が C0/C1 の定義と一致

### Audit 4: Fixed-F direct comparison ✅ PASS
- Contiguous cuts: Γ_F^exact = 4.0, C_dyn^(0.1) = 3.991, Δ = 0.009
- B2 best: Γ_F^exact = 2.616, C_dyn^(0.1) = 2.380, Δ = 0.237（τ=0.1 では期待される差）

## C3.3b 結果（再確認）

全環境で d_F ≈ 0.8:
- A2_deph_2105: avg d_F = 0.8026 ± 0.0254
- amp_damping: avg d_F = 0.8302 ± 0.0230
- mixed: avg d_F = 0.8430 ± 0.0279

## 理論的結論

**C_Γ^(0) は有限時間構造選択の予測子ではない**

これは実装エラーではなく、開量子系の物理的性質：
- Γ_F^exact(0) は瞬間的な減衰率を捉える
- 有限時間ダイナミクスは高次効果を含む
- 初期レートで予測される構造と、実際の有限時間発展で選択される構造は根本的に異なる

**C_dyn^(0) と C_dyn^(τ) は異なる時間スケールの汎関数**

Spec v2.0 の C_dyn^(0) と C_dyn^(τ) を「同じものの近似」と考えるのではなく、
独立した時間スケール汎関数として扱う必要がある。

**F* = F*(ρ, L, λ, τ) が本質的**

τ は独立パラメータであり、構造選択に本質的な役割を果たす。

## Gate 3 の再定義

### 旧 Gate 3
環境から λ を独立に決定

### 新 Gate 3
環境 (L, ρ) から (λ, τ) を独立に決定、または直接 L, ρ, τ → F*

τ の選択が構造選択に本質的な影響を与えるため、λ だけでなく τ も環境から決定する必要がある。

## Phase C の成果

1. **C0-C2.5**: open-system cost の検証と凍結
   - M₂ の open-system proxy 性を棄却
   - C_Γ^(0) = Γ_F^exact(0) を operational canonical cost として凍結
   - Spec v2.0 で一般形 Φ(F; λ, ρ, L, τ) を確立

2. **C3.1-C3.3**: λ 識別と有限時間予測の検証
   - τ_env 候補の事前定義（τ_Γ, τ_L, τ_γ）
   - matched-seed control で λ の構造選択への影響を確認
   - **C3.3b FAIL**: 初期レート予測と有限時間ダイナミクスが根本的に異なる構造を選択

3. **C3.3c**: 有限時間汎関数の監査
   - 実装の正しさを確認
   - FAIL が物理的発見であることを確定

## 次の課題

Phase D（一般性検証）の前に、以下の理論的整理が必要：

1. C_dyn^(0) と C_dyn^(τ) の理論的関係の解明
2. τ の物理的意味の明確化
3. Gate 3 の新しい定式化（λ と τ を同時に決定する方法）

## 修正されたファイル

- `analysis/gamma_F.py`: gamma_exact 関数に factor of 2 を追加
- `analysis/c3_3c_finite_time_audit.py`: 新規作成（監査スクリプト）
- `analysis_output/c3_3c_finite_time_audit.json`: 監査結果

## コミット

全 Phase C 結果をコミット準備完了。
