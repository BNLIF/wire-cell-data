# TPC Parameters Documentation

## Class Overview
`TPCParams` is a C++ class within the WCP (Wire-Cell Prototype) namespace that manages parameters and corrections for a Time Projection Chamber (TPC) detector.

## Internal Variables

### Wire Configuration Parameters
- `m_pitch_u`: Wire pitch for U plane (default: 3)
- `m_pitch_v`: Wire pitch for V plane (default: 3)
- `m_pitch_w`: Wire pitch for W plane (default: 3)
- `m_ts_width`: Time slice width (default: 3.2)
- `m_angle_u`: Angle for U plane (default: 1.0472)
- `m_angle_v`: Angle for V plane (default: -1.0472)
- `m_angle_w`: Angle for W plane (default: 0)

### Positioning Parameters
- `first_u_dis`: First U wire distance (default: 0)
- `first_v_dis`: First V wire distance (default: 0)
- `first_w_dis`: First W wire distance (default: 0)

### Configuration Parameters
- `nrebin`: Rebinning factor (default: 4)
- `time_offset`: Time offset (default: 4)
- `electron_lifetime`: Electron lifetime in ms (default: 1000)

### Correction Graphs
- `gu`, `gv`, `gw`: TGraph pointers for wire plane corrections
- `flag_corr`: Boolean for correction status

### Particle Information
- `g_proton`, `g_muon`, `g_pion`, `g_kaon`, `g_electron`: TGraph pointers for particle dQ/dx
- `g_proton_r2ke`, `g_muon_r2ke`, `g_pion_r2ke`, `g_kaon_r2ke`, `g_electron_r2ke`: TGraph pointers for range to kinetic energy

### Particle Masses (in MeV)
- `mass_proton`: 938.272 MeV
- `mass_neutron`: 939.565 MeV
- `mass_muon`: 105.658 MeV
- `mass_pion`: 139.571 MeV
- `mass_neutral_pion`: 134.977 MeV
- `mass_kaon`: 493.677 MeV
- `mass_electron`: 0.511 MeV

### Space Charge Effect (SCE) Correction
- `flag_PosEfield_corr`: Position E-field correction flag
- `h3_Dx`, `h3_Dy`, `h3_Dz`: 3D histograms for position corrections
- `h3_Ex`, `h3_Ey`, `h3_Ez`: 3D histograms for E-field corrections

## Correction and Calibration Systems

The TPCParams class implements several correction and calibration systems to account for various detector effects and ensure accurate measurements.

### Electron Lifetime Correction

#### Overview
The electron lifetime correction accounts for the attenuation of drifting electrons due to impurities in the liquid argon.

#### Implementation
```cpp
double get_attenuation_ratio(double drift_time) {
    double ratio = 1;
    if (drift_time < 0) drift_time = 0;
    if (electron_lifetime >= 1000) {
        return ratio;
    } else {
        ratio = exp(-drift_time/electron_lifetime + drift_time/200.);
        return ratio;
    }
}
```

Key Parameters:
- `electron_lifetime`: Electron lifetime in milliseconds (default: 1000 ms)
- `drift_time`: Time spent by electrons drifting to the collection plane

### Wire Plane Corrections

#### Structure
The system maintains three TGraph objects for wire plane corrections:
- `gu`: U-plane correction factors
- `gv`: V-plane correction factors
- `gw`: W-plane correction factors

#### Initialization
The `init_corr_files` method sets up wire plane corrections:
```cpp
void init_corr_files(
    TString file_u="input_data_files/calib_u_corr.txt", 
    int ndata_u = 2401,
    TString file_v="input_data_files/calib_v_corr.txt", 
    int ndata_v = 2401,
    TString file_w="input_data_files/calib_w_corr.txt", 
    int ndata_w = 3457
);
```

#### Correction Factor Calculation
The system calculates correction factors based on wire positions and angles:
```cpp
double get_corr_factor(
    WCP::Point& p,
    double offset_u, double slope_yu, double slope_zu,
    double offset_v, double slope_yv, double slope_zv,
    double offset_w, double slope_yw, double slope_zw
);
```

### Geometric Calibration

#### Wire Pitch Configuration
Manages the spacing between wires in each plane:
- U plane: `m_pitch_u` (default: 3 mm)
- V plane: `m_pitch_v` (default: 3 mm)
- W plane: `m_pitch_w` (default: 3 mm)

#### Wire Angle Configuration
Controls the orientation of wire planes:
- U plane: `m_angle_u` (default: 1.0472 radians)
- V plane: `m_angle_v` (default: -1.0472 radians)
- W plane: `m_angle_w` (default: 0 radians)

#### First Wire Positions
Maintains the starting positions of each wire plane:
- U plane: `first_u_dis`
- V plane: `first_v_dis`
- W plane: `first_w_dis`

### Time-Based Corrections

#### Time Slice Configuration
- `m_ts_width`: Width of time slices (default: 3.2)
- `time_offset`: Timing offset (default: 4)
- `nrebin`: Rebinning factor (default: 4)

#### Drift Time Corrections
Incorporates electron lifetime effects into signal processing:
1. Calculates attenuation based on drift time
2. Applies correction factors to charge measurements
3. Accounts for readout timing effects

### Combined Correction Application

#### Sequential Correction Process
1. Apply geometric corrections
   - Wire pitch adjustments
   - Angular corrections
   - Position-dependent factors

2. Time-based corrections
   - Drift time effects
   - Electron lifetime
   - Timing alignment

3. Signal corrections
   - Wire plane response
   - Cross-talk effects
   - Gain variations

#### Correction Flags and Control
- `flag_corr`: Master switch for corrections
- `flag_PosEfield_corr`: Space charge effect corrections
- Individual enable/disable options for specific corrections

### Calibration Data Management

#### File Handling
- Supports multiple input file formats (ROOT, txt)
- Handles both binary and ASCII calibration data
- Maintains separate calibration files for different correction types

#### Data Validation
- Checks for file existence and readability
- Validates calibration data ranges
- Ensures consistency across different correction types

#### Memory Management
- Proper allocation and deallocation of calibration data
- Resource cleanup in destructor
- Dynamic updating of calibration information

### Quality Control

#### Validation Checks
1. Range checking for all correction factors
2. Boundary condition enforcement
3. Error handling for invalid inputs

#### Diagnostic Output
- Prints correction factors for sample points
- Reports maximum and minimum corrections
- Provides statistical summaries of corrections

### Integration with Other Systems

#### Space Charge Effect Integration
- Combines geometric corrections with SCE
- Adjusts correction factors based on field distortions
- Maintains consistency between different correction types

#### Particle ID Integration
- Applies corrections before PID calculations
- Accounts for local field variations
- Maintains calibration for different particle types

## Particle Identification (PID) System

The particle identification system in TPCParams enables the identification and characterization of different particle types through their energy deposition patterns and ranges in the detector.

### Supported Particle Types
The system handles five primary particle types:
1. Protons (mass: 938.272 MeV)
2. Muons (mass: 105.658 MeV)
3. Pions (mass: 139.571 MeV)
   - Charged pions
   - Neutral pions (mass: 134.977 MeV)
4. Kaons (mass: 493.677 MeV)
5. Electrons (mass: 0.511 MeV)

### Data Storage Structure
The PID system maintains two sets of calibration graphs:

#### dQ/dx Calibration Graphs
TGraph objects storing the charge deposition patterns:
- `g_proton`: Proton dQ/dx profile
- `g_muon`: Muon dQ/dx profile
- `g_pion`: Pion dQ/dx profile
- `g_kaon`: Kaon dQ/dx profile
- `g_electron`: Electron dQ/dx profile

#### Range to Kinetic Energy Graphs
TGraph objects for range-energy conversions:
- `g_proton_r2ke`: Proton range to kinetic energy mapping
- `g_muon_r2ke`: Muon range to kinetic energy mapping
- `g_pion_r2ke`: Pion range to kinetic energy mapping
- `g_kaon_r2ke`: Kaon range to kinetic energy mapping
- `g_electron_r2ke`: Electron range to kinetic energy mapping

### Energy Deposition Models

#### ArgoNeut Recombination Model
The system implements the ArgoNeut recombination model for converting between dE/dx and dQ/dx:

1. dQ/dx to dE/dx Conversion:
```cpp
double func_dEdx_from_dQdx(double dQdx, double e, double alpha, double beta) {
    double result = exp(dQdx * beta * 23.6e-6/1.38/e) - alpha;
    result = result/(beta/1.38/e);
    return result;
}
```

2. dE/dx to dQ/dx Conversion:
```cpp
double func_dQdx_from_dEdx_by_ArgoNeut_model(double dEdx, double e, double alpha, double beta) {
    double result = log(alpha + beta/1.38/e*dEdx)/(23.6e-6*beta/1.38/e);
    return result;
}
```

Where:
- e: Electric field strength
- alpha: ArgoNeut model parameter (default: 0.93)
- beta: ArgoNeut model parameter (default: 0.212)
- 23.6e-6: Ionization energy in liquid argon (eV)

### Initialization and Management

#### Calibration Data Loading
The `init_PID_dq_dx` function manages the initialization of PID calibration data:
1. Cleans up any existing calibration graphs
2. Loads dQ/dx profiles from the first input file
3. Loads range-to-kinetic-energy mappings from the second input file
4. Associates each particle type with its corresponding calibration data

Example usage:
```cpp
init_PID_dq_dx("stopping_ave_dQ_dx.root", "ave_range_to_kenergy.root");
```

#### Memory Management
The class implements careful memory management for calibration data:
1. Automatic cleanup in the destructor
2. Proper deletion of existing graphs before loading new ones
3. Null pointer checks before operations

### Access Methods
The system provides getter methods for accessing calibration data:

#### dQ/dx Profile Access
```cpp
TGraph* get_pion_dq_dx();
TGraph* get_proton_dq_dx();
TGraph* get_muon_dq_dx();
TGraph* get_kaon_dq_dx();
TGraph* get_electron_dq_dx();
```

#### Range to Kinetic Energy Access
```cpp
TGraph* get_pion_r2ke();
TGraph* get_proton_r2ke();
TGraph* get_muon_r2ke();
TGraph* get_kaon_r2ke();
TGraph* get_electron_r2ke();
```

#### Mass Access
```cpp
double get_mass_proton();
double get_mass_neutron();
double get_mass_kaon();
double get_mass_pion();
double get_mass_muon();
double get_mass_pi0();
double get_mass_electron();
```

### Integration with Space Charge Effect
The PID system works in conjunction with the Space Charge Effect corrections to provide accurate particle identification in regions of distorted electric field:

1. Local electric field corrections are applied first
2. dQ/dx measurements are corrected for field variations
3. Particle identification is performed using the corrected values

## Space Charge Effect (SCE) Correction System

The Space Charge Effect is a phenomenon in Time Projection Chambers where the buildup of slow-moving positive ions creates distortions in the electric field and consequently affects both the drift paths of electrons and the local electric field magnitude.

### SCE Data Structure
The class maintains six 3D histograms to handle SCE corrections:
- Position Correction Histograms:
  - `h3_Dx`: X-direction position correction
  - `h3_Dy`: Y-direction position correction
  - `h3_Dz`: Z-direction position correction
- Electric Field Correction Histograms:
  - `h3_Ex`: X-component E-field correction
  - `h3_Ey`: Y-component E-field correction
  - `h3_Ez`: Z-component E-field correction

### Detector Dimensions
The SCE correction system uses the following detector dimensions:
- X dimension: 256 cm
- Y dimension: 232.5 cm
- Z dimension: 1037 cm

### Scaling Factors
The correction system employs specific scaling factors for each dimension:
- X scale: 2.50/2.56
- Y scale: 2.50/2.33
- Z scale: 10.0/10.37

### Correction Implementation

#### Position Correction
The `func_pos_SCE_correction` function applies position corrections by:
1. Converting real-space coordinates to correction map coordinates
2. Applying boundary conditions (0.001 to max-0.001)
3. Scaling coordinates appropriately
4. Interpolating correction values from 3D histograms
5. Applying corrections with proper scaling

#### Electric Field Correction
The E-field correction process:
1. Uses a nominal field (E0) of 0.2739 kV/cm
2. Applies the ArgoNeut recombination model with parameters:
   - α (alpha) = 0.93
   - β (beta) = 0.212
3. Corrects the electric field components:
   - Ex = E0 - correction_x
   - Ey = correction_y
   - Ez = correction_z

#### dQ/dx Correction
The `func_dQdx_after_Pos_Efield_SCE_correction` function:
1. Calculates the local electric field magnitude
2. Converts dQ/dx to dE/dx using the current field
3. Converts back to dQ/dx using the nominal field
4. Applies the ArgoNeut recombination model

### Initialization Process
The `init_Pos_Efield_SCE_correction` function:
1. Loads correction maps from an input ROOT file
2. Creates new 3D histograms with proper binning
3. Transfers correction data to the new histograms
4. Validates the transfer by comparing bin contents
5. Prints diagnostic information including:
   - Sample point values
   - Maximum/minimum correction values
   - Histogram dimensions and ranges

### Safety Features
- Boundary checking for all coordinates
- Fallback values for out-of-bounds points
- Validation of correction magnitudes
- Proper memory management for histograms

## Internal Functions

### Constructor/Destructor
- `TPCParams()`: Constructor that initializes default values
- `~TPCParams()`: Destructor that cleans up TGraph and histogram pointers

### Getters/Setters
#### Wire Configuration
- `set_pitch_u/v/w(double)`, `get_pitch_u/v/w()`: Wire pitch manipulation
- `set_angle_u/v/w(double)`, `get_angle_u/v/w()`: Wire angle manipulation
- `set_first_u/v/w_dis(double)`, `get_first_u/v/w_dis()`: First wire distance manipulation

#### Time Configuration
- `set_ts_width(double)`, `get_ts_width()`: Time slice width
- `set_time_offset(int)`, `get_time_offset()`: Time offset
- `set_nrebin(int)`, `get_nrebin()`: Rebin factor

#### Particle Mass Access
- `get_mass_proton/neutron/kaon/pion/muon/pi0/electron()`: Particle mass getters

### Correction Functions
- `init_corr_files(TString, int, TString, int, TString, int)`: Initialize correction files
- `get_corr_factor(Point&, double...)`: Calculate correction factor
- `get_attenuation_ratio(double)`: Calculate attenuation ratio based on drift time

### Space Charge Effect (SCE) Functions
- `init_Pos_Efield_SCE_correction(TString)`: Initialize SCE corrections
- `func_pos_SCE_correction(Point&)`: Apply position SCE correction
- `func_dx_after_Pos_Efield_SCE_correction(...)`: Calculate dx after SCE correction
- `func_dQdx_after_Pos_Efield_SCE_correction(...)`: Calculate dQ/dx after SCE correction

### Particle Identification Functions
- `init_PID_dq_dx(TString, TString)`: Initialize particle ID dQ/dx data
- `get_pion/proton/muon/kaon/electron_dq_dx()`: Get particle dQ/dx graphs
- `get_pion/proton/muon/kaon/electron_r2ke()`: Get particle range-to-KE graphs
- `func_dQdx_from_dEdx_by_ArgoNeut_model(...)`: Convert dE/dx to dQ/dx using ArgoNeut model
- `func_dEdx_from_dQdx(...)`: Convert dQ/dx to dE/dx