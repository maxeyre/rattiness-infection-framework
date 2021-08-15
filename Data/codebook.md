# Codebook: Linking rattiness, geography and environmental degradation to spillover *Leptospira* infections in marginalised urban settings
## Rat data
- `point_id`: unique sampling location ID
- `valley`: Pau da Lima valley ID (1/2/4)
- `X` and `Y`: coordinates (EPSG:32724 - WGS 84 / UTM zone 24S)
- `data_type`: rat index
- `outcome`: index outcome (traps: +ve/-ve, plates: number +ve, burrows: number, faeces: presence 0/1, trails: presence 0/1)
- `offset`: number of trials (traps: 1 or 0.5 if closed early, plates: number remaining, burrows: 1, faeces: 1, trails: 1)
- `offset_req`: traps: 1 if trapped closed early, 0 otherwise. Ignore for other indices.
- `rel_elev`: elevation relative to bottom of valley (m)
- `dist_trash`: 3D shortest distance to nearest refuse dump (m)
- `lc_30m_imperv`: Proportion of land cover which is impervious within a 30m radius

## Human data
- `idnova`: unique participant ID
- `ID`: unique household ID
- `outcome`: serological evidence of *Leptospira* infection (0/1)
- `X` and `Y`: coordinates (EPSG:32724 - WGS 84 / UTM zone 24S)
- `age`: current age in years
- `gender`: gender (0 female, 1 male)
- `valley`: Pau da Lima valley ID (1/2/4)
- `income_pcap`: household income per capita (USD$/day)
- `educ_yrs`: years of full-time education
- `lc_30m_imperv`: Proportion of land cover which is impervious within a 30m radius
- `elev_levels`: household elevation level (1: low 0-6.7m, 2: medium (6.7-15.6m), 3: high (>15.6m))
- `hillside`: household located on a hillside
- `rel_elev`: elevation relative to bottom of valley (m)
- `open_sewer_10m`: open sewer within 10m of household
- `exposed_to_sewer`: no barrier between household and sewer
- `work_constr`: work in construction
- `work_salesman`: work as travelling salesperson
- `work_trashcollect`: work in refuse collection
- `work_mud`: work involves contact with mud
- `work_floodwater`: work involves contact with floodwater
- `work_sewer`: work involves contact with sewer water
- `floodwater_freq_bin`: contact with floodwater in last 6 months (0: never/rarely, 1: sometimes, 2: frequently)
- `sew_water_freq_bin`: contact with sewer water in last 6 months (0: never/rarely, 1: sometimes, 2: frequently)
- `pred.mean.rattiness`: mean predicted rattiness at household location
- `illiteracy`: functionally illiterate adult (0: no, 1: yes)
