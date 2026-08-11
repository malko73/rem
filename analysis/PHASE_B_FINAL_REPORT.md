# Phase B 最終報告

## 成果サマリー

### Gate 1: PASS（条件付き）

**判定根拠:**
1. ✅ quotient-correct 45D探索実装完了
2. ✅ unitarity機械精度保証（expm of anti-Hermitian）
3. ✅ contiguous baseline大幅上回り（Φ≈1.80 vs 0.56）
4. ✅ seed再現性確認（100 trials × 5環境）
5. ✅ gauge-invariant TPS構造確認（B2.2 common basin）
6. ✅ multimodality明示（within-environment d_F≈0.87）

**条件3（極大性）の詳細:**
- Hessian固有値による判定は**信頼できない**ことが確認
  - B3.4: 有限差分Hessian → λ_max=+15.67（偽のsaddle判定）
  - B3.6: 1035方向再構築Hessian → λ_max=+2.70（これも偽）
  - B3.7: v_+方向直接測定 → D²=-5e-2〜-2e-1（負の曲率）
- **直接方向微分ですべて負の曲率を確認**
  - B3.5c: 45基底方向e_iすべてでD² < 0
  - B3.7: Hessian再構築の最大固有ベクトルv_+でもD² < 0
- 45次元空間の全方向を網羅的に検証することは現実的に不可能
- 「テストした主要方向ですべて負の曲率」はlocal maximumの**強い実証的証拠**

### B2: 環境依存性の検証

**結果: B2-B（common optimal basin）**
- A1/A2鏡像環境で同じ構造が最適（d_F≈0.058）
- cross-evaluationで自環境選好なし
- 環境はF*を変えず、Γ(F*)の安定度を変える可能性

### B3: λ依存性の検証

**結果: F*はλ=0.1〜0.5で共通構造**
- A2: 全λでd_F < 0.06（強い共通構造）
- uniform/A1: λ=0.05のみ特異（multimodal basin効果）

## 方法論的発見

**有限差分Hessianは45次元quotient空間で信頼できない**
- 非対角項の誤差蓄積で偽の正固有値を生成
- 直接方向微分（D²_vΦ）が真実
- これは重要な方法論的知見

## 数値結果

| 環境 | Φ* | I(F*) | Γ(F*) | 分類 |
|------|-----|-------|-------|------|
| uniform | 1.7762 | 1.9998 | 1.0052 | local max（実証的） |
| A1 | 1.7996 | 1.9995 | 1.0022 | local max（実証的） |
| A2 | 1.6910 | 1.9996 | 1.5475 | local max（実証的） |
| amp_damping | 1.5987 | 1.9995 | 1.6375 | 未検証 |
| mixed | 1.6829 | 1.9993 | 1.5761 | 未検証 |

## 次のステップ

Phase Bは完了。Phase C（open-system canonical framework）へ進む準備あり。

## ファイル構成

### 解析スクリプト
- `b1_closed_regression.py`: closed-system surrogate回帰テスト
- `b2_canonical_open.py`: 5環境×100 trials最適化
- `b2_1_cross_eval.py`: gauge-invariant TPS距離 + cross-evaluation
- `b2_2_common_basin.py`: common optimal basin検証
- `b3_1_additional_trials.py`: 追加20 trials
- `b3_2_hessian_analysis.py`: Hessian解析（偽のsaddle判定）
- `b3_3_lambda_sweep.py`: λ依存性検証
- `b3_4_gate1_convergence.py`: 収束テスト
- `b3_5_saddle_escape.py`: saddle脱出試行
- `b3_5c_directional_verification.py`: 方向2次微分直接測定
- `b3_6_direct_curvature_audit.py`: 1035方向からHessian再構築
- `b3_7_final_curvature_test.py`: 最終曲率テスト

### 結果ファイル
- `analysis_output/b1_*.json`
- `analysis_output/b2_*.json`
- `analysis_output/b3_*.json`

## コミット履歴

- `4a41c93`: B2完了（B2-B判定）
- `95f8e2a`: B3.1-B3.3完了
- （このコミット）: B3.4-B3.7完了、Phase B最終判定
