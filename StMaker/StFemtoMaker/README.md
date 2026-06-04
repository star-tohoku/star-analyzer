# StFemtoMaker naming rules

Canonical location for femto **species keys**, **particle keys**, and **channel names**.
When adding particles or channels, follow these rules and update `FemtoConfig` YAML accordingly.

## Species keys (`species` map in femto YAML)

- Format: `lower_snake_case`
- One key per candidate source registered in the event store
- Examples: `proton`, `anti_proton`, `kaon_plus`, `phi`
- Each species entry requires:
  - `builderType`: `track` or `resonance`
  - `particleKey`: builder dispatch key (see below)

## Particle keys (`particleKey` field)

- Selects builder logic inside generic `TrackPidBuilder` / `ResonanceBuilder`
- Phase 1 supported values: `proton` (track), `phi` (resonance from KK)
- Use PDG-inspired names; charge sign goes in species key when needed (`anti_proton`)

## Channel names

- Format: `{partA}_{partB}` using species keys
- Example: `phi_proton` with `partA: phi`, `partB: proton`
- Histogram suffixes use channel name: `hKstarSE_phi_proton`, `hCF_phi_proton`

## YAML flat-key layout

`YamlParser` is flat; femto config uses prefixed keys:

```yaml
speciesKeys: proton,phi
species_proton_builderType: track
species_proton_particleKey: proton
species_phi_builderType: resonance
species_phi_particleKey: phi

channelName: phi_proton
partA: phi
partB: proton
signalMin: 1.01
signalMax: 1.03
normQMin: 0.2
normQMax: 0.3
```

## Agent / developer checklist

1. Add species block with unique species key
2. Implement or extend builder dispatch for new `particleKey` in `StFemtoMaker`
3. Add channel entry referencing existing species keys
4. Add histogram entries to `config/hist/hist_anaFemto*.yaml` using channel name suffix
5. Document new keys in this file

See also `docs/femto_track_sharing_concerns.md` for overlap risks when extending species.
