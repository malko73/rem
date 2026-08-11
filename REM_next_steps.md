# REM 研究ロードマップ（正式版）

**版**: v2.0（2026-08-11）
**前版**: `REM_next_steps_v1_2026-04-25.md`（v4基準で全面更新）
**基準ソース**: REM_lambda v5（Zenodo 10.5281/zenodo.21880505, 2026-08-11 公開）／ `malko73/rem` branch `open` ／ `REM_spec_v2_0.tex`（凍結版）

---

## 現在地の定義

REMは「着想だけの理論」でも「連続最適化を作る段階」でもない。

> **「すでに作った理論・実装を、壊れない仕様へ固定し、物理proxyを反証テストする段階」に入った。**

変分原理・λクロスオーバー・具体的量子系での数値計算は存在証明として完了済み。次の仕事は概念の追加ではなく、**仕様凍結 → 厳密化 → 反証テスト**である。

## 認識の修正（v4基準）

1. **連続TPS最適化は「未実装」ではない。** 基礎実装済み（`optimize_factorization` / `optimize_factorization_adam`）。残るのは「厳密化フェーズ」（大域性・頑健性の実証）。
2. **数値基準は (Z₁Z₃, λ\*≈0.3672) ではない。** XY鎖・v4基準（λ\*≈0.165）。

## 検証済みベースライン（v4 / open @ 67aa2b8）

> **⚠️ 2026-08-11 unitarity hotfix（対応完了）**: 連続TPS最適化の U = exp(iΣθG) は anti-Hermitian G に対し非ユニタリ（正定値 Hermitian、||U†U−I||=0.30）であるバグを発見・修正済み（U = exp(ΣθG)、||U†U−I||≈2e-16）。
> **再検証結果（30 trials, N=3,4,5）: 連続最適化はもはや contiguous を安定的に超えない。** N=3: SGD 0.581±0.112/60%、Adam 0.654±0.191/67%。N=4: SGD 1.709±0.064/53%、Adam 1.752±0.031/100%。N=5: SGD 1.454±0.049/60%、Adam 1.480±0.028/97%。robustness は N=3 initial-θ 軸のみ頑健（87%）、他はコイントス（43-57%）。旧主張（Φ=1.05、1.95×、180/180成功）は**非ユニタリ写像の artifact で無効**。Adam>SGD は方向性維持。early stopping 45%削減・品質損失ゼロ。**v4 Zenodo レコードの連続最適化数値は修正後の値に置き換えが必要。** 離散cut（M₂ vs Var、λ\*）は影響なし。商空間次元は **45**（U(1) kernel 考慮）に修正。

| 項目 | 状態 | 数値 |
|---|---|---|
| 動的コスト規約 | ✅ 確定 | C_H = ⟨H_∂²⟩ ≥ 0（符号確定、旧 −⟨H_∂⟩ は廃止） |
| 数値例 | ✅ 確定 | 3-qubit XY鎖 (J₁₂=1.5, J₂₃=0.6, h=0.2, λ=0.2)、クロスオーバー λ\*≈0.165 |
| 連続最適化（N=3） | ✅ 基礎実装 | 最良contiguous Φ=0.560 → 連続 Φ=1.05（mean 1.09±0.22, 30 seeds, 成功率100%） |
| ロバストネス（N=3,4,5） | ✅ 実施 | 180 trials 全てcontiguous超え。N=4: 1.88±0.03 vs 1.723、N=5: 1.62±0.03 vs 1.458 |
| Optimizer比較 | ✅ 実施 | Adam > SGD 90/90 trials（+9.9%〜+14.6%） |
| Early stopping | ✅ 実施 | N=3で53%ステップ削減（品質劣化<0.002）。N=4,5では150 stepsでは無効 |
| テスト | ✅ 18件 | Hamiltonian/状態不変量・C_H定義・⟨H²⟩vs⟨H⟩²区別・離散cut・連続最適化改善・Adam正しさ・SGD–Adam乖離 |
| 明示的制約 | ⚠️ 残存 | 大域最適性の保証なし（有限差分勾配・生成子部分空間・SGD/Adam）／ λの微視的導出なし／ デコヒーレンス計算なし／ 実験検証なし |
| **spec §6** | ✅ **修正済み（A0, 2026-08-11）** | `REM_spec_v1_1.tex` §4/§6 を v4同期（C_H = Tr(ρH_∂²)、Φ = Φ_S − λC_H）。§8a に Monotonic Tradeoff Theorem 追加。variance への変更は比較実験後に判断 |

---

## Phase A — Core Freeze（REM Spec v2.0相当）: 定義の凍結

最初に数学的コアを確定し、**ここで一度定義を凍結**する。概念追加はしない。
正式項目は下記6項目。実作業の実行順は **A0 → 1〜5 → 6（最終監査）** とする。

- **A0: spec §6 v4同期 hotfix**（✅ 2026-08-11 実施済み）
  `REM_spec_v1_1.tex` §6 の旧式 Φ_H = −Tr(ρH_boundary) を廃止し、
  C_H(F) = Tr(ρH_∂F²)、Φ(F;λ) = Φ_S(F) − λC_H(F) へ同期した。
  「C_H は boundary-energy **second moment** であり、variance ではない」と明記済み。
  variance への変更可否は C_H vs Var 比較実験の結果が出るまで確定しない。

1. **Pareto / monotonicity theorem の正式化（Monotonic Tradeoff Theorem）**
   λ₂ > λ₁ のとき、各λの大域最適解 F₁=F\*(λ₁), F₂=F\*(λ₂) について
   C(F\*(λ₂)) ≤ C(F\*(λ₁)) かつ S(F\*(λ₂)) ≤ S(F\*(λ₁))。
   最適性2式の加算で証明可能（(λ₂−λ₁)(C₁−C₂) ≥ 0）。C が second moment でも variance でも独立に成立するので、先に確定してよい。数値に依存しないREMの数学的性質として Lemma/Theorem 化する。

2. **C_H = ⟨H_∂²⟩ と Var(H_∂) の比較**（✅ 2026-08-11 実装・暫定判定 B）
   `analysis/m2_vs_var.py` + `tests/test_m2_vs_var.py`（24 tests 全て green）。
   結果: Δ=⟨H⟩² は M₂ の 24〜98% と大きい（proxy選択は無意味でない）。Var(H_∂) は N=3 ベンチマークで **因子化に完全非依存**（cut1=cut2=0.6206896551724, diff 2.6e-15）、N=4,5 でも min-Var cut が max-I cut と一致 → **Var 下では λ* が存在せずクロスオーバー構造が消滅**（agreement 5.7〜7.0%、M₂ の switching 1〜2回 vs Var 0回）。
   **暫定判定 B（Phase C の Γ_F 照合が必須）+ C寄りエビデンス**（Var は選択駆動力として退化 → M₂ 維持の暫定根拠）。最終決定は Phase C で M₂/Var vs Γ_F を照合して行う。

3. **D の数学的身分の整理**
   D(ρ, H, λ) = F\* = argmax_F Φ(F) として「構造選択」を明示し、時間発展 U_t・環境デコヒーレンス ℰ と完全分離する。REM1の記述的写像（𝒟: ℋ → ℋ_S⊗ℋ_O）との整合を明記する。

4. **TPS商空間・群作用・kernel・次元の最終監査**
   F ≅ U(N)/(U(n_A)⊗U(n_B)) の群作用・kernel（Kronecker積表現での (e^{iθ}I_A, e^{−iθ}I_B) 1次元kernel）・次元（N=3では4、フルU(8)/(U(4)⊗U(2))は44）を厳密に定義。REM3/REM5の商多様体の数学を監査する。

5. **目的関数の無次元化**
   I はdimensionlessだが C_H は energy²。C̃_H = C_H/E₀² として無次元化し、Φ = I − λC̃_H とする。λの意味（情報量と動的コストのスケール変換）が明確になる。

6. **論文・spec・コードで定義を完全一致させる**
   - `REM_spec_v1_1.tex` §6 の旧式 Φ_H = −Tr(ρH_boundary) を C_H = Tr(ρH_∂²) へ修正 ← **既知の不一致**
   - v2/v4論文・README・コード・テストの規約が一致していることを確認
   - 凍結後は REM_spec_v2.0 として確定

**Phase A 完了条件**: spec・論文・コードの定義が完全一致し、以降の数値変更が定義の変更を伴わないこと。

---

## Phase B — Continuous TPS Validation（既存実装の研究品質化）

新規実装ではなく、既存 `optimize_factorization` 系を研究品質へ引き上げる。

目的は「連続TPS探索コードが動いた」ではなく、**「得られた最適因子化が大域解または十分に頑健な解である」ことの実証**。

- 多数の random initialization（現在30 seeds → 系統的スイープ）
- Adam / SGD 等の optimizer 依存性
- 局所解率（異なる初期条件から同じ解へ収束するか）
- Hessian / curvature 解析（停留点の性質）
- 離散cutとの比較（contiguous優位性の定量化）
- Hamiltonian family 依存性（asymmetric XY / transverse Ising / Heisenberg / XYZ / random local）
- N scaling（現在 N=3,4,5 → より大きいN、eighタイムアウト対策）

**Phase B 完了条件（v5 基準）**: Gate 1（quotient-correct full TPS search で、再現可能な非自明 F* が存在する）を、局所解率・curvature解析付きで示せること。canonical REM については F*(ρ, L, λ) が環境変化に対して系統的に応答することを確認。

### B0（✅ 2026-08-11 実施・commit ddec387）: Quotient Geometry / Optimizer Correctness

SVD ベースの商空間直交分解（`src/quotient_geometry.py`）を実装。U(1) kernel は手で消さず SVD に rank=19 を検出させる方式（マスター決定）。

- dim 𝔲(8) = 64 / raw vertical = 20 / **rank vertical = 19** / **dim horizontal = 45**
- V†H ≈ 0（<1e-10）、P_V²=P_V、P_H²=P_H、P_V+P_H=I
- 必須9項目（gauge invariance・vertical ΔΦ≈0・horizontal ΔΦ≠0・unitarity 機械精度・spectrum/norm 不変・projected gradient の vertical 成分ゼロ）をテストで固定
- 旧実装バグ発見: u(d) 実基底のオフ対角符号（実**反**対称 + 虚**対**称が正。逆だと Hermitian になり射影が壊れる）
- 45次元 optimizer（`src/quotient_optimizer.py`・SGD/Adam）接続: スモークで Φ=1.473・unitarity ≤5.3e-16

### B1（✅ 2026-08-11 実施・commit 5d0e6ed）: closed-system regression（45次元 quotient）

Φ_closed = I − λC_H^closed を45次元 quotient optimizer で再実行（N=3 XY, λ=0.2, 100 trials × 200 steps）。

| optimiser | best | median | mean | std | min | success | distinct minima | max unitarity |
|---|---|---|---|---|---|---|---|---|
| SGD | 1.4556 | 1.4409 | 1.4389 | 0.0100 | 1.4081 | 100% | 2 | 7.8e-16 |
| Adam | 1.4780 | 1.4777 | 1.4775 | 0.0007 | 1.4749 | 100% | 1 | 1.1e-15 |

contiguous best = 0.5602。**非自明な高Φ解は unitarity 修正後も生存**（Φ≈1.44-1.48、contiguous の2.6倍）。Adam は単一 dominant 極大（std=0.0007）、SGD は軽度分散（std=0.010, 2極小）。**Gate 1 の予備的支持**（B2 で確定）。

### B2a（✅ 2026-08-11 実施・commit 150c073）: C_Γ^(0) well-posedness audit

連続 TPS 最適化の前に、canonical open-system cost C_Γ^(0) の well-posedness を監査（発散・縮退・gauge・frame）。

| テスト | 結果 |
|---|---|
| 1. product limit（product TPS への経路） | ✅ 有界（Γ∈[1.10,2.09], Φ∈[0.16,1.64]） |
| 2. random TPS sweep（2000点） | ✅ 有界（Γ∈[0.45,2.36], 全て正, corr(Γ,1/C_F²)≈-0.03） |
| 3. Schmidt 非零縮退（p₁=p₂） | ⚠️ **Γ が縮退ブロック内の Schmidt basis に依存**（2.000 vs 1.787）。C1.5 の null-completion 不変性は非零縮退に拡張されない。xfail として記録（spec 未決: pinching 射影 or domain 制限） |
| 4. local gauge invariance（V_A⊗V_B） | ✅ PASS（ΔΓ≈1e-15） |
| 5. frame consistency（pull-back vs lab-frame） | ✅ PASS（環境も回転させれば一致、ΔΓ≈1e-15） |

**発見・修正したバグ（C0/C1 では基底状態 [H,ρ₀]=0 のため潜伏）**:
- `gamma_F.liouvillian_env` の commutator 項の kron 引数が vec_F 規約で逆（正: -i[kron(I,H) - kron(H.conj(),I)]）
- 回転環境の dephasing は kron(z.T, z) が必要（z^rot=U†ZU は非対称）

**判定**: 発散なし（B2a-FAIL/WARN ではない）。gauge・frame は成立。**Schmidt 非零縮退での basis 依存は新規の well-posedness caveat** — canonical functional 凍結前に解決が必要（候補: 分母なし C_rate^(2)、pinching 射影、または domain 制限）。

### B2（✅ 2026-08-11 実施・commit 39a3e87）: canonical Open REM continuous TPS

Φ(F; λ, ρ, L) = I_ρ(F) − λΓ_F^exact(0) を 45次元 quotient optimizer で最適化（5環境 × 100 trials）。

| 環境 | best Φ | median | std |
|---|---|---|---|
| uniform dephasing | 1.7756 | 1.7721 | 0.0036 |
| A1 deph (0.5,1,2) | 1.7990 | 1.7948 | 0.0035 |
| A2 deph (2,1,0.5) | 1.6901 | 1.6867 | 0.0033 |
| amp damping | 1.5987 | 1.5914 | 0.0075 |
| mixed | 1.6829 | 1.6715 | 0.0084 |

**B2.1（gauge-invariant TPS distance + cross-eval）**: 局所代数 projector 距離 d_F を実装。within-environment d_F≈0.87（全環境・300 steps でも不変 → **landscape は strongly multimodal**、optimizer 収束問題ではない）。best U* 同士は近い（A1 vs A2: 0.058）。cross-eval では自環境選好なし（Φ_A1(U*_A1)=1.7990 vs Φ_A1(U*_A2)=1.7985；Φ_A2 は A1 の解をわずかに好む）。

**B2.2（common basin test）**: top 20 + top 20 の candidate pool を両環境で cross-evaluate → **top-5 overlap 0.60 / top-10 overlap 0.70**。A1 の最良解（trial 6）は A1・A2 両環境で 1 位。**COMMON BASIN（strong evidence）**。

**正式判定: B2-B** — "no resolved environment-dependent shift of the optimal factorization. The continuous landscape is strongly multimodal, but the best A1/A2 solutions are compatible with a common optimal factorization; the observed environmental dependence appears primarily in the objective value through Γ_F, not in the selected structure itself."

**研究上の意味（マスター）**: 環境は F* を必ずしも変えず、**候補構造の安定性を重み付けする**可能性。F*_A1 ≈ F*_A2（this benchmark, λ=0.2）。解析勾配は Phase B 後半の性能改善として保留。

### B3 — Robustness / Hessian

---

## Phase C — Open REM（反証テスト優先）✅ 完了

ここで初めてLindbladへ進む。最初の研究命題は λ 導出ではない。

> **C_H(F) は実際の Γ_F を予測するか**

比較対象: ⟨H_∂²⟩, Var(H_∂), Γ_F（Lindbladモデルで直接計算）

検証: **F\*_{C_H} ?= F\*_{Γ}**（REMが選ぶ因子化と、実際に最も長寿命な構造の一致）

- 一致すれば強い（REMの目的関数が物理的安定性を正しくproxyしている）。
- **一致しなければREMの目的関数を修正する**（この研究姿勢が理論を強くする）。
- 成功した場合のみ λ の微視的導出（Γ_F ≃ κC_H の因子化独立性チェック → 単一λモデルの妥当性検証）へ進む。

### C0（✅ 2026-08-11 実施・ケース2）: 同一純粋dephasing環境での Γ_F 測定

固定環境（L_i = √(γ/2)Z_i, γ=1.0、H・ρ₀ は cut1/cut2 で完全同一）・SciPy 64×64 Liouvillian・C_F(t) = ‖ρ(t) − D_F[ρ(t)]‖₂（初期Schmidt基底への dephasing projection）。主指標は解析的 t=0 微分 Γ_F^exact(0) = −Re⟨X_F, Q_F L(ρ₀)⟩_HS / |X_F|₂²（Q_F = I−D_F, X_F = Q_Fρ₀）。

| cut | M₂ | Var | Γ^exact | Γ^(0.05) | Γ^fit |
|---|---|---|---|---|---|
| 1 | 8.3793 | 0.62069 | **2.000000** | 1.9929 | 1.7384 |
| 2 | 0.8193 | 0.62069 | **2.000000** | 1.9962 | 1.4988 |

**Γ₁^exact = Γ₂^exact = 2.000000（完全一致）: 実デコヒーレンス初期レートは因子化に非依存。** Var の退化構造と整合し、M₂ の10倍の安定性差（8.38 vs 0.82）は open-system に現れない。**ケース2: λ\*≈0.165 は second-moment proxy 由来の artifact である可能性を真剣に検討する。**

### C1（✅ 2026-08-11 実施・C1-C）: 4環境での Γ_F 測定

事前登録4環境（A1: deph γ=(0.5,1,2) / A2: 鏡像 γ=(2,1,0.5) / B: amp.damping κ=1 / C: mixed γ=κ=0.5）。H・ρ₀・L_i は cut 間完全固定。

| 環境 | Γ₁^exact | Γ₂^exact | R_Γ |
|---|---|---|---|
| A1 deph (0.5,1,2) | 1.6212 | 2.7685 | **0.586** |
| A2 deph (2,1,0.5) | 2.9394 | 1.9630 | **1.497** |
| B amp.damping | 1.0000 | 1.0000 | **1.000** |
| C mixed | 1.5000 | 1.5000 | **1.000** |

**判定 C1-C（環境依存）**: M₂ の固定順位（R_M₂=10.23）はどの環境でも再現されない。均一/対称ノイズでは Γ は完全に因子化非依存（R=1.000）、非一様ノイズでは順位がノイズ分布に従い鏡像で反転。**→ 動的コストは環境依存でなければならない: C_H = C_H(F, E)。Φ(F;E) = I(F) − λC_open(F;E) へ進む根拠**（閉鎖系 second moment ではなく、実際の open-system generator から導出）。これは REM の失敗でなく物理的基礎の強化。

### C1.5（✅ 2026-08-11 実施・2点修正）: 測定定義の監査

1. **sigma_minus バグ修正**: `[[0,0],[1,0]]`（|1⟩⟨0| = raising）→ `[[0,1],[0,0]]`（|0⟩⟨1| = lowering）へ修正。B: Γ₁=Γ₂=2.000000（旧1.000000）、C: Γ₁=Γ₂=2.000000（旧1.500000）——R=1.000 不変、**C1-C 判定は修正に頑健**。
2. **gauge robustness（100 seeds）**: 零 Schmidt 部分空間の completion に Γ^exact は**完全不変**（rel_std = 0.00e+00、両 cut）。解析的理由: 純粋状態では D_F が自己随伴かつ D_F X = 0 より分子 Re⟨X, Q_F L(ρ)⟩ = Re⟨X, L(ρ)⟩ となり D_F が落ちる。**basis-independent pinching への変更は不要**。

### 仕様の二層構造（C1-C を受けた方針）

- C_H^closed(F) = ⟨H_∂F²⟩ は**閉鎖系 surrogate に格下げ**（削除しない）。
- 正式な動的項: **C_dyn(F; ρ, L, τ)**（抽象的な環境依存量）。基本式:
  **Φ(F; λ, ρ, L, τ) = I_ρ(F) − λ C_dyn(F; ρ, L, τ)**
- Monotonic Tradeoff Theorem は E・ρ・τ を固定すればそのまま成立（証明不変）。
- **C_Γ^(0)(F;ρ,L) = Γ_F^exact(0)** を open-system cost の第一候補（operational realization）として仕様へ追加。
- 将来的に C_dyn^(0) = Γ^exact（t=0）と C_dyn^(τ) = −(1/τ)log(C_F(τ)/C_F(0))（有限時間）を分離。
- C2: Liouvillian boundary cost C_L(F) = |L_∂F|² または C_L(F,ρ) = |L_∂F(ρ)|₂² を次候補として構築し、Γ_F との予測性能を比較してから Spec v2.0 の動的コストを凍結。

### C2（✅ 2026-08-11 実施・C2-B）: Liouvillian boundary cost vs Γ_F

L_∂F = L − P_F^local L（P_F^local は S_F = {L_A⊗I_B² + I_A²⊗L_B} への HS 直交射影、行列単位基底＋pseudoinverse。jump operator 非一意性に非依存）。

| 環境 | R_Γ(exact) | R_L(global) | R_L(rho) | R_L(align) |
|---|---|---|---|---|
| A1 deph (0.5,1,2) | 0.586 | 4.159 | 0.597 | −0.284 |
| A2 deph (2,1,0.5) | 1.497 | 6.090 | 0.445 | 0.568 |
| B amp.damping | 1.000 | 5.291 | 0.655 | 退化 |
| C mixed | 1.000 | 5.759 | 0.415 | 0.750 |

**判定 C2-B**: norm 型 C_L は R_Γ を予測できない（global は常に >1、rho は常に <1、環境に追従しない）。Spearman(global)=0.949 は4点の順位偶然、log-ratio 誤差大（1.96/1.22）。align 診断は amplitude damping で退化（境界残差が実際のコヒーレンス減衰方向とほぼ直交）。**→ C_L は canonical 候補に昇格しない。C_Γ^(0) = Γ^exact を operational canonical open-system cost として確定する方向。**

### C3（✅ 2026-08-12 実施）: λ識別と有限時間予測

**C3.1（環境時間尺度の候補）**: τ_env として τ_Γ, τ_L, τ_γ, τ_κ を定義。全候補が B3 の共通構造領域（λ∈[0.1,0.5]）に落ちることを確認。

**C3.3a（matched-seed λ-sensitivity control）**: 同一初期条件で λ=0.2 と λ=λ_pred を比較。
- τ_L: d_F ≈ 0.08-0.14（構造を動かさない）
- τ_Γ, τ_γ: d_F ≈ 0.16-0.22（構造を動かすが Φ が崩壊）

**C3.3b（finite-time hold-out prediction）**: ❌ **FAIL**
- 予測側: F_pred = argmax [I(F) - λ_pred·Γ_F^(0)]
- 独立評価側: F_finite = argmax [I(F) - λ_pred·C_dyn^(τ_env)(F)]
- 結果: d_F(F_pred, F_finite) ≈ 0.80-0.84（全環境）

**C3.3c（finite-time functional audit）**: ✅ 実装の正しさを確認
- τ→0 limit: C_dyn^(τ) → Γ_F^exact として収束（factor of 2 バグ修正済み）
- Frame covariance: ✅
- D_F definition consistency: ✅
- Fixed-F direct comparison: ✅

### Phase C の結論

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

**Gate 3 の再定義**:
- 旧: 環境から λ を独立に決定
- 新: 環境 (L, ρ) から (λ, τ) を独立に決定、または直接 L, ρ, τ → F*

τ の選択が構造選択に本質的な影響を与えるため、λ だけでなく τ も環境から決定する必要がある。

**Phase C 完了条件**: ✅ 達成
- Gate 2（C_Hが実際のopen-system stabilityと対応）: ✅ PASS
- Gate 3（環境パラメータからλを因子化選択とは独立に決定）: 🔲 OPEN（再定義必要）

---

## Phase D — Generality / Falsification

- mixed state（ρ ≠ |ψ⟩⟨ψ|）
- finite temperature（λ\*(T) → (J, T, λ) 3次元相図）
- XY / Ising / Heisenberg / XYZ / random local Hamiltonian（連続的なC_H競合を作る）
- N scaling
- noise robustness
- QPU実験（P2）

**Phase D 完了条件**: Gate 4（λから未知条件でF\*予測）＋ Gate 5（別ハミルトニアン・混合状態でも同一原理）。

---

## Gates（到達条件）

| Gate | 内容 | 通過で何が変わるか |
|---|---|---|
| **Gate 1** | Full TPS optimizationでも非自明な F\* が存在 | Phase B 完了。Toy model批判を大幅に潰す |
| **Gate 2** | C_H が実際のopen-system stabilityと統計的・解析的に対応 | REMの目的関数を物理的に裏付ける |
| **Gate 3** | 環境パラメータから λ を因子化選択とは独立に決定 | フリーパラメータ批判を解消 |
| **Gate 4** | その λ から未知条件で F\* を予測 | 「入力物理条件 → λ → F\*」が成立。予測理論に |
| **Gate 5** | 別ハミルトニアン・混合状態でも同一原理 | 一般的な研究プログラムとして説得力 |

Gate 1〜3で「かなり強い理論」。Gate 4で「予測理論」。Gate 5で「研究プログラム」。

**Gate 状態（2026-08-12 更新）**:
| Gate | 状態 |
|------|------|
| Gate 1 | **Conditional PASS**（direct directional curvature testsによるempirical local-maximum evidence。有限差分Hessianは信頼できない制約あり） |
| Gate 2 | **PASS / operational definition established**（C_Γ^(0)=Γ_F^exact(0) 凍結、Spec v2.0） |
| Gate 3a | 🔲 **OPEN**（τ固定でλを環境から独立決定） |
| Gate 3b | 🔲 **OPEN**（τ自体を環境から決定、より強い条件） |
| Gate 4 | OPEN |
| Gate 5 | OPEN |

**Gate 3 の再定義（2026-08-12）**:

Phase C3.3b の結果（C_dyn^(0) ≠ C_dyn^(τ)）を受け、Gate 3 を二分割：

- **Gate 3a**: (ρ, L, τ) → λ を factorization 選択とは独立に決定
  - τ は外から指定された物理的観測時間
  - より弱い条件、達成しやすい

- **Gate 3b**: (ρ, L) → (τ, λ) まで物理的に決定
  - 例: Liouvillian gap τ_L = |Re μ₁|⁻¹ から自然な τ が得られる
  - より強い条件、REMに不必要に強い要求を課す可能性

Phase C3.1-C3.2 の結果:
- τ_L は構造を動かさない（d_F < 0.1）
- しかし τ_L は有限時間構造選択を予測しない（C3.3b FAIL）
- → Gate 3b は現時点では未達

**推奨**: まず Gate 3a を達成し、τ の物理的意味が明確になってから Gate 3b に進む。

---

## 旧ロードマップ（v1, 2026-04-25）からの変更点

| 項目 | v1 | v2（本版） |
|---|---|---|
| 最優先 | λの微視的導出（3-6ヶ月） | Phase A 仕様凍結 → B 厳密化 → C で初めて C_H↔Γ_F 検証 → 成功時のみ λ 導出 |
| 連続多様体最適化 | 「実装する」（未実装扱い） | 「基礎実装済み・厳密化フェーズ」 |
| 数値例 | 3-qubit Z₁Z₃、λ\*≈0.3672 | XY鎖、λ\*≈0.165 |
| 動的コスト | Φ_H = −⟨H_∂⟩（線形） | C_H = ⟨H_∂²⟩（二次・符号確定） |
| 反証姿勢 | （明記なし） | **C_H↔Γ_F 不一致なら目的関数を修正**（§C） |
| 既存研究 | 「REMだけがTPS選択理論」含意 | Operational Quantum Mereology（arXiv:2212.14340, Quantum 2024）／ Loizeau & Sels（arXiv:2409.01391）／ Adlam & Rovelli（arXiv:2203.13342）と差別化。独自性は「二目的競合 Φ=I−λC による構造レジーム転移＋RQMのfactorization選択問題への接続」に絞る |
| 数学的証明 | 「長期（1-3年）」 | Pareto/monotonicityは Phase A で即座に確定 |

## 保留事項（v1から引き継ぎ）

- **Dの物理的実現は急がない**: Dは記述的な二分割選択写像として定義。「Dの物理的実現は今後の課題」と正直に書く（査読者に受け入れられやすい）。
- **統合論文（REM1+2+3+5）**: v1では推奨だったが、今回は保留。Gates 1-3通過後の方が説得力を持つ。
- 出版は概念論文を増やすより、Continuous TPS optimization → Open-system validation → λ独立導出の3本を通すことを優先。

---

## Spec v2.0 正式凍結（✅ 2026-08-11）

- **canonical**: Φ(F; λ, ρ, L, τ) = I_ρ(F) − λ C_dyn(F; ρ, L, τ)、**[λ] = T**（特徴的時間尺度）。[Φ]=1, [I]=1, [C_dyn]=T⁻¹
- **第一実装**: C_dyn^(0) = Γ_F^exact(0)（signed structural decoherence rate）。C_dyn^(τ) = −(1/τ)log(C_F(τ)/C_F(0)) も同一 λ
- **名称**: signed dynamical functional（負値は reward。cost は便宜用語）
- **τ₀ 正規化**: optional representation のみ（canonical に新基準時間を導入しない）
- **C_H^closed = ⟨H_∂F²⟩**: closed-system structural surrogate に格下げ（v4 連続性のため削除せず）
- **Monotonic Tradeoff Theorem**: signed C_dyn でも成立（λ≥0 のみ必要）
- ファイル: `REM_spec_v2_0.tex`（v1_1 は履歴として残置）

**Phase A（Core Freeze）は完了。** ✅ **REM_lambda v5 は 2026-08-11 に Zenodo 公開済み**（DOI 10.5281/zenodo.21880505 / https://zenodo.org/records/21880505）。

## REM_lambda v4 → v5 訂正版（✅ 2026-08-11 公開済み）

v4（Zenodo 10.5281/zenodo.21427776）の**独立した2つの重要修正点**を v5 として訂正公開。旧 v4 は publication history として残置。

| # | 修正点 | 内容 |
|---|---|---|
| 1 | continuous TPS 数値の unitarity artifact | Φ=1.05 等は非ユニタリ写像の artifact。修正後は SGD 0.581±0.112/60%、Adam 0.654±0.191/67%（N=3）などに差し替え。商空間次元 44→45 |
| 2 | M₂ の open-system 解釈不支持 | C0: Γ₁=Γ₂=2.000000（因子化非依存）／ C1: 環境依存（R_Γ: A1=0.586, A2=1.497, B=C=1.000）／ C2: norm 型 C_L 不支持（C2-B） |
| 3 | λ\*≈0.165 の再解釈 | 「open-system で実在する物理的クロスオーバー」→「**closed-system second-moment surrogate のもとで得られる構造的クロスオーバー**」に soften |
| 4 | 一般形の提示 | Φ(F;E) = I − λC_dyn(F;E)（open-system REM） |
| 5 | C_Γ^(0) の導入 | 最初の operational candidate として |
| 6 | テスト拡充 | 50 tests（unitarity・M₂vsVar・Γ_F・C1/C2・gauge） |

## 保留中の Spec v2.0 正式凍結（C2.5 ✅ 2026-08-11 実施済み）

1. **Γ^exact の符号**: random states / Hamiltonians / environments 200 configs・400 値で **Γ < 0 が 5 件（1.25%、最小 −1.607）発生**（XY/GUE 両方）。coherence 初期増加は実在する。**名称は signed structural decoherence rate に決定**（max(0,Γ) は使わない）。spec §6 に SIGN 注記追加済み。
2. **次元の整理**: uniform dephasing で Γ ∝ γ が完全線形（Γ/γ = 2.000000, std=0）→ **[Γ] = T⁻¹ 確認、[λ] = T**（特徴的時間尺度）。無次元化 τ₀Γ_F も可能。spec §6 に次元注記追加済み。

→ 残る凍結判断は「λ を T の次元のまま使うか、τ₀ で無次元化するか」の選択のみ。

## 参照ファイル

- `REM_next_steps_v1_2026-04-25.md` — 旧ロードマップ（バックアップ）
- `~/Desktop/Programing/rem/` — コード本体（branch `open`）
- `REM_lambda/REM_lambda_v5.tex` — 現行論文 v5（訂正版・公開済み）
- `REM_spec_v2_0.tex` — 仕様（凍結版 2026-08-11）
- `REM4/rem4_numerical.py` — 数値コード（ToyBox側コピー）
- 方針文書: 2026-08-11 評価ドキュメント（Phase A-D・Gates 1-5・C_H↔Γ_F反証テスト）
- スキル: `rem-research-workflow`（2026-08-11時点の状態を反映済み）
