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

굽힘 반경 **R** 과 측면 오프셋 **dy** 만 주면, **동일한 두 개의 원호**로
S-bend(offset bend)를 만들어 **끝단 진행각이 항상 시작각(0°)으로 복귀**하는
세그먼트입니다.

- 첫 번째 호: 곡률 `+s/R`, `+s·θ` 회전
- 두 번째 호: 곡률 `−s/R`, `−s·θ` 회전 → 시작 heading 으로 복귀
- `s = sign(dy)`

즉 시작 포트와 동일한 방향을 보면서 측면으로 `dy` 만큼 평행 이동한
끝 포트를 자동 생성합니다. SBEND 와 달리 곡률 반경이 고정값으로
직접 지정되므로, 최소 굽힘 반경 제약(곡률 손실)을 보장하기 좋습니다.

## 물리 / 기하 공식

```
dy = 2 R (1 − cos θ)        →     θ = arccos(1 − |dy| / (2R))
전진거리  dx = 2 R sin θ      (R, dy 로 자동 결정)
```

- 각 호의 회전각 **θ** 는 R 과 dy 로 완전히 결정됩니다.
- 도달 가능 범위: `|dy| ≤ 4R` (θ ≤ 180°). 초과 입력 시 `|dy| = 4R` 로
  클램프되어 θ = 180° 가 됩니다.
- `dy = 0` 또는 `R ≤ 0` 이면 straight stub(퇴화) 으로 처리됩니다.

## 파라미터

| key | 라벨 | 기본값 | 의미 |
|---|---|---|---|
| `radius` | Bend radius [µm] | 500.0 | 두 원호의 굽힘 반경 R |
| `dy` | Lateral offset dy [µm] (±) | 100.0 | 목표 측면 오프셋 (부호 = 좌/우) |

(`width` / `xoffset` / `yoffset` / `angleoffset` 등 공통 오프셋 파라미터도
다른 세그먼트와 동일하게 적용됨)

## Auto-Offset 와의 관계

`rel_bent` 는 시작·끝 heading 이 동일(0° 복귀)하므로 `sbend` 와 동일하게
**straight-like** 로 취급됩니다(`_segment_signed_curvature` 가 양 끝단 κ=0
반환). 따라서 rev25 의 auto bend-offset 접합 보정에서 추가 offset 을
유발하지 않습니다.

## 검증

- 끝점 y 가 모든 도달 가능 케이스에서 목표 `dy` 와 정확히 일치
  (예: R=500, dy=+100 → end=(435.89, 100.00), θ=25.84°).
- 클램프 케이스(R=100, dy=500 > 4R=400) → dy=400, θ=180° 로 정상 제한.
- 끝 heading 은 항상 시작각(`end_a = a0`).
- `python3 -m py_compile` 구문 검사 통과.

## 하위 호환

rev26 이전 프로젝트 파일(JSON)에는 `rel_bent` 세그먼트가 없으므로
영향이 없으며, 기존 세그먼트 종류·동작은 완전히 보존됩니다.
