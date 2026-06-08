# rev27 — REL_BENT (offset bend) 세그먼트 추가 요약

## 변경점 (rev26 → rev27)

| 위치 | 추가/수정 |
|---|---|
| 헤더 docstring | 버전 rev27, 세그먼트 목록에 `REL_BENT` 추가, rev27 기능 소개 |
| `Segment` dataclass 주석 | kind 목록에 `'rel_bent'` 추가 |
| `segment_geometry()` | `kind == "rel_bent"` geometry 블록 신설 (두 원호 합성) |
| `SEG_PARAMS["rel_bent"]` | `radius`, `dy` 두 파라미터 등록 (Add Segment 폼 노출) |
| `SEG_DISPLAY["rel_bent"]` | 드롭다운 표시명 추가 |
| `_segment_signed_curvature()` | straight-like 그룹에 `rel_bent` 추가 (양 끝 κ=0 취급) |
| `_show_about()` | About 다이얼로그에 rev27 / REL_BENT 표기 |

## REL_BENT 란?

굽힘 반경 **R** 과 **절대 y 타깃 `dy`** 만 주면, **끝단 진행각이 항상
절대 0°(+x 방향)** 가 되도록 하면서 끝점이 정확히 `y = dy` 에 도달하는
오프셋 bend 세그먼트입니다. **들어오는 각도와 무관**하게 동작합니다.

3개의 반경-R 원호로 분해 구성됩니다:

1. **헤딩 교정 호**: 단일 원호로 들어오는 heading `a0` 를 `0°` 로 펴줌
   (이때 생기는 수직 이동 = `Δy_A`)
2. **대칭 S-bend (2호)**: 남은 수직 변위 `H_res = (dy − y0) − Δy_A` 를
   채우는 두 원호 (heading 0 → ±θ → 0)

a0 = 0(수평 진입)이면 1단계는 생략되어 단순 대칭 S-bend 가 됩니다.
SBEND 와 달리 곡률 반경이 고정값으로 지정되어 최소 굽힘 반경(곡률 손실)
제약을 보장하기 좋습니다.

## 절대 y 타깃 + `yoffset(n)`

`dy` 는 **상대 변위가 아니라 도달할 절대 y 좌표**입니다. 따라서
`dy = yoffset(5)` 라고 쓰면 **segment 5 의 끝점 절대 y** 에 정확히
정렬되는 bend 가 만들어집니다.

`yoffset(n)` 은 모든 숫자 입력 필드/파라미터 식에서 쓸 수 있는 헬퍼로,
segment n 의 끝점 포트 절대 y 를 돌려줍니다(참조 대상 segment 가
시퀀스상 먼저 있어야 함). `Project.eval_globals()` 로 평가 namespace 에
주입되며, geometry 파라미터를 해석하는 모든 지점(rebuild/redraw/add/
edit/GDS export/hit-test/Measure)에서 동작합니다.

## 물리 / 기하 공식

```
H      = dy − y0                                  (필요한 전체 수직 이동)
Δy_A   = R·sign(−a0)·(cos a0 − 1)                 (헤딩 교정 호의 수직 이동)
H_res  = H − Δy_A                                 (S-bend 가 채울 잔여 변위)
θ      = arccos(1 − |H_res| / (2R))               (S-bend 각 호 각도)
```

- 헤딩 프로파일(global)을 `a0 → 0 → ±θ → 0` 으로 적분하여 중심선 생성,
  끝 heading 은 **항상 정확히 0°**.
- 도달 가능 범위: `|H_res| ≤ 4R`. 초과 시 `4R` 로 클램프됩니다.
- a0 ≈ 0 이고 dy == y0 이면 퇴화(길이 0) 처리됩니다.

## 파라미터

| key | 라벨 | 기본값 | 의미 |
|---|---|---|---|
| `radius` | Bend radius [µm] | 500.0 | 각 원호의 굽힘 반경 R |
| `dy` | Target y [µm] (abs; e.g. yoffset(5)) | 100.0 | 도달할 **절대 y 좌표** |

(`width` / `xoffset` / `yoffset` / `angleoffset` 등 공통 오프셋 파라미터도
다른 세그먼트와 동일하게 적용됨. 단, rel_bent 의 `dy` 는 절대 y 타깃이므로
공통 `yoffset` 파라미터로 시작점을 옮기면 그만큼 H 가 재계산됨)

## Auto-Offset 와의 관계

`rel_bent` 는 `sbend` 와 동일하게 **straight-like** 로 취급됩니다
(`_segment_signed_curvature` 가 양 끝단 κ=0
반환). 따라서 rev25 의 auto bend-offset 접합 보정에서 추가 offset 을
유발하지 않습니다.

## 검증

- **끝 heading 이 들어오는 각도(a0=0/30/−25/40/90/45°)와 무관하게 항상
  정확히 절대 0°**, 동시에 끝점 y 가 목표 절대 y `dy` 와 정확히 일치.
  (예: a0=40°, y0=20, dy=20 → end_y=20.000, head=0° / a0=90°, dy=300 →
  end_y=300.000, head=0°).
- `yoffset(n)`: segment n 끝점 y 반환, 식 `dy="yoffset(5)"` 로 정렬 동작,
  미존재 segment 는 명확한 에러.
- 클램프 케이스(|H_res| > 4R) → 4R 로 정상 제한.
- `python3 -m py_compile` 구문 검사 통과.

## 하위 호환

rev26 이전 프로젝트 파일(JSON)에는 `rel_bent` 세그먼트가 없으므로
영향이 없으며, 기존 세그먼트 종류·동작은 완전히 보존됩니다.
