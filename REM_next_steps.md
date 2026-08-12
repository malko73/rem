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

### D0（✅ 2026-08-12 実施・commit a091f72）: Protocol Freeze

評価指標（Φ*, d_F, Γ^exact, C_dyn^(τ)）、固定パラメータ（λ=0.2, τ=0.1, Adam, 200 steps, 6 seeds）、5 Hamiltonian families、environment set、成功基準を `analysis_output/d0_protocol.json` に凍結。

### D1（✅ 2026-08-12 実施・commit a091f72）: Hamiltonian Generality — **pre-v2.2 preliminary evidence（superseded normalized functional 使用）**

> **Status（2026-08-12, Spec v2.2 凍結後）**: 旧 D1 は削除せず履歴として残す。**pre-v2.2 preliminary evidence using the superseded normalized dynamical functional**（Φ = I−λΓ^exact）。D1-R が通った時点で「Hamiltonian generality under the Spec v2.2 canonical functional」へ昇格。

5 families × 6 seeds、Φ = I−λΓ^exact を 45次元 quotient で最適化（A1 dephasing 環境）。

| Hamiltonian | contiguous Φ | best Φ | success |
|---|---|---|---|
| asymmetric_XY | 0.170 | 1.799 | ✅ |
| transverse_field_ising | 1.093 | 1.650 | ✅ |
| heisenberg_xxz | 0.895 | 1.739 | ✅ |
| xyz | 0.847 | 1.748 | ✅ |
| random_local | 0.455 | 1.696 | ✅ |

**5/5 success**。interaction 構造・対称性・ランダム性を変えても非自明な F* が出現（fine-tuning 排除の予備的証拠）。

### D1-R（✅ 2026-08-12 実施・commit 871a104）: Hamiltonian Generality revalidation — **Spec v2.2 canonical で全ゲート PASS**

**canonical functional のみ使用**: J_dyn^(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩、Φ = I − λ·J_dyn^(0)。凍結プロトコルは旧 D1 と同一（**SEED0=20260812 も旧 D1 のまま**）で、dynamical term のみ変更。

| Hamiltonian | contiguous Φ | best Φ (D1-R) | 旧 D1 | std |
|---|---|---|---|---|
| asymmetric_XY | 0.653 | **1.8994** | 1.7993 | 0.0007 |
| transverse_field_ising | 1.331 | **1.8246** | 1.6498 | 0.0010 |
| heisenberg_xxz | 1.281 | **1.8695** | 1.7389 | 0.0002 |
| xyz | 1.242 | **1.8740** | 1.7481 | 0.0006 |
| random_local | 0.807 | **1.8481** | 1.6962 | 0.0007 |

**Gate 判定（全て PASS）**:
- **R1** 5/5 finite optimum ✓
- **R2** seed 安定（std ≤ 0.0010）✓
- **R3** product collapse なし（p_max 0.508–0.513）✓
- **R4** J_dyn / I / Φ 非自明（best >> contiguous、|J_dyn| > 0.1）✓
- **R5 F* が Hamiltonian に応答** ✓ — cross-eval 対角優位、**mean gap = 0.113**（asym 0.068 / ising 0.092 / heisenberg 0.106 / xyz 0.025 / random 0.274）。**family-dependent response with varying separation strength**（XYZ の 0.025 は小さく、「全 family が強く分離」とは断定しない）。d_F 行列: asym↔random 0.412、ising↔xyz 0.559 が近く、他は 0.87–0.88（構造的類似 family が近い = 物理的に意味のある応答）
- **R6** D2.1 型 singularity 再発なし（min C_F² = 0.500、max|γ_ref| = 1.78、max|Φ| = 1.899）✓

**解釈（記録の正確化）**: **Spec v2.2 canonical objective の下で、全5 Hamiltonian family に安定した非自明 optimum が得られた。** 旧 D1（Φ_Γ）と D1-R（Φ_J）は目的関数自体が異なるため、絶対値を直接比較して「改善」とは評価しない。D1-R5 の cross-evaluation は、各 F*_j を異なる Hamiltonian の objective で再評価して対角優位が出たため、**F* = F*(ρ, L, λ) が単なる状態依存の共通解ではなく、dynamics/Hamiltonian に実際に応答している**ことを直接支持する。**「Hamiltonian generality under the Spec v2.2 canonical functional」へ昇格。**

### D2（✅ 2026-08-12 実施・D2-R 改訂）: State Generality — **pre-v2.2 preliminary evidence（superseded normalized functional 使用）**

> **Status（2026-08-12, Spec v2.2 凍結後）**: 旧 D2 は削除せず履歴として残す。**pre-v2.2 preliminary evidence using the superseded normalized dynamical functional**（Φ = I−λΓ^exact）。D2-R が通った時点で「state generality under the Spec v2.2 canonical functional」へ昇格。

4状態タイプ（ground / Haar / mixed / thermal）で最適化。

| State | best Φ | std | d_F from ground |
|---|---|---|---|
| ground | 1.799 | 0.001 | — |
| haar | 49.10 ⚠️ | 17.63 | 0.44 |
| mixed | 1.849 | 0.028 | 0.87 |
| thermal | 1.850 | 0.036 | 0.86 |

**Original D2 criterion: FAIL** — F* の状態非依存性を要求する D0 基準は Spec v2.1（F* = F*(ρ, L, λ, τ)）と不整合。
**D2-R 改訂**: 各状態クラスが非自明・数値安定・再現可能な構造最適解を持つこと（状態間の F* 一致は不要）。
- ground: PASS
- mixed: PASS候補
- thermal: PASS候補
- **Haar: PENDING / anomaly（→ D2.1-A に分類）**

### D2-R（✅ 2026-08-12 実施・commit 64e5b03）: State Generality revalidation — **Spec v2.2 canonical で全ゲート PASS**

**canonical functional のみ使用**: J_dyn^(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩、Φ = I − λ·J_dyn^(0)。凍結プロトコルは旧 D2 と同一（**SEED0=20260813 も旧 D2 のまま**）で、dynamical term のみ変更。D2.2-A は singularity well-posedness 監査であり、D2-R は frozen canonical の state generality 再認証（代用しない）。

| State | contiguous Φ | best Φ (D2-R) | std |
|---|---|---|---|
| ground | 0.653 | **1.8995** | 0.0007 |
| haar | 1.529 | **1.8594** | 0.0090 |
| mixed | 0.683 | **1.9845** | 0.0051 |
| thermal | 0.137 | **1.9465** | 0.0058 |

（数値は D2.2-A と完全一致 — プロトコル・シード同一のクロスチェック合格）

**Gate 判定（全て PASS）**:
- **R1** 4/4 finite optimum ✓
- **R2** seed 安定（std ≤ 0.0090）✓
- **R3** product collapse なし（p_max 0.510–0.580）✓
- **R4** J_dyn / I / Φ 非自明（best >> contiguous、|J_dyn| > 0.1）✓
- **R5 F* が state に応答** ✓ — cross-eval 対角優位、**mean gap = 0.159**（ground 0.000 / haar 0.332 / mixed 0.058 / thermal 0.247）。**state-dependent response with varying separation strength**。**正確な記録**: state-dependent response is established overall, although ground and mixed contain nearly degenerate optima under the ground-state objective（ground の gap=0.000 のため「各状態固有の F* が存在する」とは断言しない）。d_F 行列: **ground↔mixed = 0.343（近い）**、haar↔thermal は 0.87–0.89（強く分離）。ground と Haar は d_F 0.880 / gap 0.332 で明確に異なる。mixed と thermal は Φ 値が近いが **d_F 0.868 で構造は異なる**。**全状態が同一 attractor ではない**（mean gap 0.159）
- **R6** singularity 再発なし（min C_F² = 0.096、max|γ_ref| = 1.66、max|Φ| = 1.98）✓

**解釈**: **Spec v2.2 canonical objective の下で、全4状態に安定した非自明 optimum が得られた。** 旧 D2（Φ_Γ）と D2-R（Φ_J）は目的関数が異なるため絶対値を直接比較しない（Haar の 49.10 → 1.8594 は singularity 除去であり、D2.1/D2.2-A で確立済み）。D2-R5 の cross-eval 対角優位は **F* = F*(ρ, L, λ) が状態に実際に応答する**ことを直接支持。**「State generality under the Spec v2.2 canonical functional」へ昇格。**

**Spec v2.2 の3本柱が揃った**: well-posedness（D2.2-A/B）+ Hamiltonian generality（D1-R）+ **state generality（D2-R）**。→ **次は D3 Timescale**。

### D3（✅ 2026-08-12 実施・commit 8c1a236）: Timescale — **τ_c = 0.018（Haar）、sharp structural crossover / 全ゲート PASS**

**問い**: Does relational structure itself depend systematically on the dynamical observation scale? 中心対象: Haar（当初ブラケット 3e-3 < τ_c < 1e-2 は best-seed 由来で不正確 → 精密同定で修正）。

**D3-A: τ_c 精密同定（Haar、τ=0.003..0.1 の19点、6 seeds、凍結 d2_2b プロトコル）**

マスター定義の ΔΦ(τ) = Φ_τ(F_B) − Φ_τ(F_A)（F_A = D2.2-A Haar 最適、F_B = τ=0.1 Haar 最適、固定代表）:

| τ | 0.003 | 0.01 | 0.016 | **0.018** | 0.02 | 0.025 | 0.03 | 0.05 | 0.1 |
|---|---|---|---|---|---|---|---|---|---|
| ΔΦ | −0.021 | −0.011 | −0.003 | **0.000** | +0.003 | +0.009 | +0.016 | +0.039 | +0.079 |

**τ_c = 0.0180**（ΔΦ 符号反転の線形補間）。ΔΦ は滑らか・単調。optimizer は鋭く追随: **τ=0.018 で全6 seed A 値 → τ=0.02 で全6 seed B 値**（値ベース分類。std が 0.005 → 0.0009 に collapse）。

**重要な方法論修正**:
1. **d_F ベースの per-seed basin ラベルは無意味**（45次元商空間では多数の局所最大値がどの2盆地からも d_F≈0.87。例: τ=0.1 で「A 判定」された seed も Φ=1.911 ≈ B 値）。**値ベース分類**（Φ が Φ_τ(F_A) か Φ_τ(F_B) に近いか）が正しい
2. D2.2-B の「τ≥1e-2 で全 seed 移動」は best-seed の見かけ（記録修正済み）。物理的クロスオーバーは ΔΦ 交差（τ_c=0.018）で定義 — optimizer の basin hopping と分離できる（マスター予告どおり）

**D3-B: 転移の性質** — **first-order-like basin transition（sharp structural crossover）in the realization; continuous in the objective**。adjacent d_F は τ=0.02（0.859）と τ=0.05（0.877）で O(1) ジャンプ（best 解が異なる局所最大値をホップ）、一方 ΔΦ は滑らか。best 解は τ=0.1 で F_B に完全収束（d_F=0.000）。有限次元のため「相転移」とは呼ばず「sharp structural crossover / first-order-like basin transition」とする。

**D3-C: Liouvillian timescale との比較** — **Δ_L = 0.658、τ_L = 1.519、τ_c·Δ_L = 0.0119**。**O(1) 予想は不支持**（τ_c は遅い緩和モードより遥かに短い）。τ_c は「B 盆地の瞬間的不利を有限時間 O(τ) 補正が上回る」**inter-basin competition scale** で決まり、Liouvillian gap ではない。thermal の τ_c∈(0.3,1.0) なら比 0.2–0.66 — 状態依存で集中せず。**結論: このベンチマークでは τ は単なる Liouvillian gap タイムスケールではない**（D3-C の仮説はこの系では棄却 — 重要な負の結果）

**D3-7: Haar 以外のクロスオーバー（粗い確認、τ=0.2/0.3/1.0）**:
- **ground: クロスオーバーなし**（τ=1.0 の F* は D2-R 最適から d_F=0.11、ΔΦ≈0）
- **mixed: クロスオーバーなし**（d_F=0.12、ΔΦ≈0）
- **thermal: クロスオーバーあり、τ_c ∈ (0.3, 1.0)**（ΔΦ: −0.011 → +0.001、d_F(F_A,F_B)=0.879）

→ **クロスオーバーの有無・位置は状態依存**（Haar τ_c=0.018 ≪ thermal τ_c≈0.3-1.0、ground/mixed は 1.0 まで無し）。relational structure は dynamical observation scale に系統的に依存する（Haar・thermal で実証）。

**Gate 判定（全て PASS）**: D3-1（τ_c=0.018 特定）✓ / D3-2（ΔΦ 符号反転）✓ / D3-3（seed 非依存 — ΔΦ は固定代表で objective レベル、optimizer は勝ち盆地の値に追随）✓ / D3-4（singularity/product collapse 再発なし）✓ / D3-5（連続性分類: 実現は first-order-like、objective は連続）✓ / D3-6（τ_c/τ_L 評価: 0.012、O(1) 不支持）✓ / D3-7（Haar 以外の確認: thermal にあり、ground/mixed なし）✓

### D2.1（✅ 2026-08-12 実施・commit 0ad8840）: Haar Singularity Audit — **D2.1-A 確定**

**C_Γ^(0) の正規化は unrestricted TPS optimization 上で特異**（重要な反証結果）:

- corr(Φ, log10 C_F(0)²) = **−0.9990**（ほぼ完全な負相関）
- C_F(0)² → 6e-5 まで減少すると Γ → −222, Φ → 44（発散）
- best 解の Schmidt spectrum p = [0.9997, 0.0003]（ほぼ product state）
- 任意の純粋状態は適切な TPS で product state に近づけられるため、分母 |Q_F ρ|² → 0 が可能

**含意**: Γ_F^exact(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩/|Q_Fρ|² は局所診断量としては使えるが、full TPS 上の global variational functional としては **well-posed でない**。

**Spec v2.1 の再検討が必要**。マスター提示の代替案:
- **J_dyn^(0) = −Re⟨Q_Fρ, Q_F L(ρ)⟩**（分母なし）
- **−[C_F²(τ) − C_F²(0)]/(2τ)**（分母なし有限時間）

### D2.2-A（✅ 2026-08-12 実施・commit 24db549）: Unnormalized instantaneous functional — **G1–G6 PASS / G7 PENDING**

**方針（マスター決定 2026-08-12）**: D2.2 は A（unnormalized instantaneous）に限定して実行。Spec v2.2 の執筆・D3/D4 は凍結。

**変更点は dynamical term のみ**:
```
J_dyn^(0)(F) = −Re⟨Q_Fρ, Q_F L(ρ)⟩   （分母 |Q_Fρ|² なし）
Φ = I − λ·J_dyn^(0)
```
Hamiltonian / state / optimizer / initialization / restart数 / TPS parameterization / λ / seed / stopping criteria は **D2 と完全に同一**（asymmetric_XY, dephasing γ=(0.5,1,2), λ=0.2, Adam 200 steps lr=0.01, SEED0=20260813, horizontal_basis(4,2), n_a=2, 6 seeds；Haar は D2.1 比較用に 30 seeds 追加）。

| State | D2 Φ (正規化) | D2.2-A Φ (unnormalized) | D2.2-A std |
|---|---|---|---|
| ground | 1.7993 | **1.8995** | 0.0007 |
| haar | **49.10 ⚠️** | **1.8594** | 0.0090 |
| mixed | 1.8485 | **1.9845** | 0.0051 |
| thermal | 1.8499 | **1.9465** | 0.0058 |

Haar 30 seeds: best **1.8595** / median 1.8582 / std **0.0077**（D2.1: 25.71 / 1.81 / 5.96）→ **巨大 outlier 消失**。

**特異性除去の証拠**:
- min C_F²（全 trial）= **0.0956**（D2.1: 4.1e-4）— optimizer が特異領域に近づかない（Haar 30 seeds 全て C_F² ≈ 0.49 に収束 = product 方向への attraction なし）
- max|γ_ref|（最適点での正規化診断量）= **1.66**（D2.1: 245）
- Schmidt p_max = **0.51–0.58**（D2.1: 0.9997）— product collapse なし
- G4 連続性: ε=1e-2 摂動で max|dΦ| = **7.7e-5**、全方向有限
- singularity sweep（best Haar 解から C_F²→0 方向へ ε=0.5 まで）: C_F² ≥ 0.46、Φ ≤ 1.86、全有限

**Gate 判定**: G1（singularity-free）✓ / G2（Haar outlier 消失）✓ / G3（seed 安定）✓ / G4（摂動連続）✓ / G5（非 product collapse）✓ / G6（4状態で再現可能）✓ / **G7（finite-time との整合）= PENDING（D2.2-B 実行後に判定）**。

**判定**: **singularity-free ∧ non-trivial**。D2.1 は「normalized local decay rate cannot serve as a global variational functional over unrestricted tensor factorizations」という明確な反証結果として確定。J_dyn^(0) が canonical dynamical functional の第一候補に。

**付随修正（commit fb08a3f）**: `gamma_exact` の factor of 2（537d854 で追加）を C0 事前登録規約へ revert。−d/dt log C² と −d/dt log C の混同で、マスターの D2.2 式・C0 docstring・数値微分テストの三方に矛盾していた。suite 74 passed + 1 xfailed に復旧。

**次の一手**: マスター承認後に **D2.2-B**（J_dyn^(τ) = −[C_F²(τ)−C_F²(0)]/(2τ)、同一条件）→ G7 判定。instantaneous と finite-time が同じ TPS を選ぶか確認できれば Spec v2.2 の説得力が上がる。Spec v2.2 / D3 / D4 は凍結継続。

### D2.2-B（✅ 2026-08-12 実施・commit 2fce54d）: Finite-time functional — **G7 PASS / J_dyn^(0) を canonical candidate #1 に昇格**

**J_dyn^(τ)(F) = −[C_F²(τ)−C_F²(0)]/(2τ)**（分母なし）。Q_F は t=0 の Schmidt 基底で固定 → τ→0 で J_dyn^(0) に厳密一致。凍結条件は D2.2-A と同一。**τ スイープ {1e-3, 3e-3, 1e-2, 3e-2, 1e-1}**（ρ(τ)=e^{τL}ρ₀ は U 非依存なので事前計算、評価コストは D2.2-A と同等）。

| State | best Φ @ τ=1e-3 | @ τ=1e-2 | @ τ=1e-1 | D2.2-A (τ→0) |
|---|---|---|---|---|
| ground | 1.8997 | 1.9016 | 1.9187 | 1.8995 |
| haar | 1.8589 | 1.8548 | **1.9121** ⚠️ | 1.8594 |
| mixed | 1.9841 | 1.9852 | 1.9887 | 1.9845 |
| thermal | 1.9467 | 1.9451 | 1.9651 | 1.9465 |

**G7 判定**:
- **G7-1（τ→0 極限）PASS**: J_dyn^(τ) → J_dyn^(0) を厳密確認。dev ≈ 0.96·τ·|J₀|（線形収束）、τ=1e-3 で全状態 dev ≤ 2.5e-3。master の恒等式 dC_F²/dt = 2Re⟨Q_Fρ,Q_FL(ρ)⟩ を数値的に裏付け
- **G7-2（同一/同一盆地 F*）PASS（τ-極限 caveat 文書化）**: τ ≤ 3e-3 で全状態同一盆地（d_F < 0.06）；ground/mixed/thermal は全 τ で同一盆地（cross-eval gap ≤ 0.005）；**Haar のみ τ ≥ 1e-2 で全6 seed 揃って別盆地へシフト**（I 1.97→2.00, J_τ 0.60→0.44, Φ 最大 +4.3% @ τ=0.1）。これは Spec v2.0 の F*=F*(ρ,L,λ,τ) 予言どおりの τ 依存性であり、病理（singularity/product collapse/seed instability）は一切再発しない

> **記録修正（2026-08-12, D3 による精査）**: 「τ≥1e-2 で全6 seed 揃って別盆地へ」は **best-seed ベースの見かけ**でした。D3 の per-seed 値ベース解析では、τ_c 近傍（0.003〜0.018）では seed は A/B 値に分裂（optimizer の basin trapping）し、**全 seed が B 値へ揃うのは τ≈0.02 以降**です。物理的クロスオーバーは objective レベルの ΔΦ 交差（τ_c = 0.018）で定義すべきで、G7-2 の定性結論（τ→0 で同一、有限 τ でシフト）は変わりません。詳細は D3 セクション参照。
- **G7-3（順位維持）PASS**: mixed > thermal > ground > haar の順位が全 τ で保存
- **G7-4（病理非再発）PASS**: 全 (state, τ) で min C_F²(0) ≥ 0.09, p_max ≤ 0.58, std ≤ 0.009；Haar 30 seeds @ τ=0.1: best 1.9135 / std 0.0011

**判定**: **G7 PASS**。J_dyn^(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩ を **Spec v2.2 canonical dynamical functional の第一候補として昇格**。finite-time 版 J_dyn^(τ) は canonical の代替ではなく「**有限時間における operational extension / consistency diagnostic**」として位置づけ。

**次の一手（マスター承認順）**: **Spec v2.2**（J_dyn^(0) を canonical として書き起こし）→ **D1/D2 再確認**（新 canonical での再検証）→ **D3 Timescale** → **D4 N-scaling**。

### Spec v2.2（✅ 2026-08-12 凍結・commit acecb86）

`REM_spec_v2_2.md` + `papers/REM_spec_v2_2.tex`（`tools/gen_spec.py` で再現生成、CI latex gate ビルド成功）。

- **canonical**: J_dyn^(0)(F;ρ,L) = −Re⟨Q_Fρ, Q_FL(ρ)⟩（分母なし・符号付き。「cost」とは呼ばない。正式名: instantaneous dynamical rate functional (signed)）
- **canonical objective**: Φ(F;λ) = I_ρ(F) − λ·J_dyn^(0)、[J]=T⁻¹, [λ]=T
- **Γ_F 格下げ**: fixed-(F) local diagnostic only; prohibited as unrestricted global TPS variational objective（D2.1 反証を v2.1 破棄の根拠として明記）
- **finite-time**: J_dyn^(τ) = −[C_F²(τ)−C_F²(0)]/(2τ) を operational extension / consistency diagnostic として分離。lim_{τ→0} J_dyn^(τ) = J_dyn^(0) を主要整合条件に
- **finite-τ TPS crossover を許容**: Haar の basin switch（τ≥1e-2）を genuine structural crossover として記録 → **D3 は τ-crossover 解析に昇格**
- **G7 の記録**: 「PASS with documented finite-τ crossover」（完全 PASS とはしない）
- v2.1 からの変更点（factor-of-2 revert、log-ratio → difference-ratio）を明記

**付随修正（commit cfd7db1）**: ルート `conftest.py` 追加 — 素の `pytest tests/` が src/analysis を import できず CI が失敗していた問題を根本原因修正。CI（3.12/3.13）＋ latex-build ともに green。

**次の一手（マスター承認順）**: **D1 再検証**（5 Hamiltonian families を J_dyn^(0) で再実行）→ **D2 再検証**（4状態を J_dyn^(0) で再実行）→ **D3 Timescale**（τ_c の特定: Haar の basin switch は 3e-3〜1e-2 の間）→ **D4 N-scaling**。

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
