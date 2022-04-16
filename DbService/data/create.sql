--
-- should run as admin_role, which should own all objects
SET ROLE admin_role;

--
--
CREATE SCHEMA val;
GRANT USAGE ON SCHEMA val TO PUBLIC;

CREATE SCHEMA tst;
GRANT USAGE ON SCHEMA tst TO PUBLIC;

CREATE SCHEMA trk;
GRANT USAGE ON SCHEMA trk TO PUBLIC;

CREATE SCHEMA sim;
GRANT USAGE ON SCHEMA sim TO PUBLIC;

-- rows are calibration tables
CREATE TABLE val.tables 
  (tid SERIAL, 
  name TEXT NOT NULL, 
  dbname TEXT NOT NULL, 
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT tables_pk PRIMARY KEY (tid),
  CONSTRAINT tables_unique_name   UNIQUE (name), 
  CONSTRAINT tables_unique_dbname UNIQUE (dbname)  );
GRANT SELECT ON val.tables TO PUBLIC;
GRANT INSERT ON val.tables TO manager_role;
GRANT UPDATE ON val.tables_tid_seq TO manager_role;

-- one row for each cid, a commit of a calibration table
CREATE TABLE val.calibrations
  (cid SERIAL, 
  tid INTEGER NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT calibrations_pk PRIMARY KEY (cid), 
  CONSTRAINT calibrations_tid_fk FOREIGN KEY (tid) REFERENCES val.tables(tid) );
GRANT SELECT ON val.calibrations TO PUBLIC;
GRANT INSERT ON val.calibrations TO val_role;
GRANT UPDATE ON val.calibrations_cid_seq TO val_role;

-- one row for an interval of validity
CREATE TABLE val.iovs
  (iid SERIAL, 
  cid INTEGER NOT NULL,
  start_run INT NOT NULL, 
  start_subrun INT NOT NULL, 
  end_run INT NOT NULL,
  end_subrun INT NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT iovs_pk PRIMARY KEY (iid), 
  CONSTRAINT iovs_cid_fk FOREIGN KEY (cid) REFERENCES val.calibrations(cid) );
GRANT SELECT ON val.iovs TO PUBLIC;
GRANT INSERT ON val.iovs TO val_role;
GRANT UPDATE ON val.iovs_iid_seq TO val_role;

-- one row for each group of IOVs
CREATE TABLE val.groups
  (gid SERIAL, 
   create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
   create_user TEXT NOT NULL,  
   CONSTRAINT groups_pk PRIMARY KEY (gid) );
GRANT SELECT ON val.groups TO PUBLIC;
GRANT INSERT ON val.groups TO val_role;
GRANT UPDATE ON val.groups_gid_seq TO val_role;

-- relationship between a group and iovs
CREATE TABLE val.grouplists
  (gid INTEGER, 
  iid INTEGER,
  CONSTRAINT grouplists_pk PRIMARY KEY (gid, iid), 
  CONSTRAINT grouplists_gid_fk FOREIGN KEY (gid) REFERENCES val.groups(gid), 
  CONSTRAINT grouplists_iid_fk FOREIGN KEY (iid) REFERENCES val.iovs(iid) );
GRANT SELECT ON val.grouplists TO PUBLIC;
GRANT INSERT ON val.grouplists TO val_role;

-- purposes, like PRODUCTION, CALIBRATION, TEST
CREATE TABLE val.purposes
  (pid SERIAL, 
  name TEXT NOT NULL, 
  comment TEXT NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT purposes_pk PRIMARY KEY (pid),
  CONSTRAINT purposes_unique_name UNIQUE (name)  );
GRANT SELECT ON val.purposes TO PUBLIC;
GRANT INSERT ON val.purposes TO manager_role;
GRANT UPDATE ON val.purposes_pid_seq TO manager_role;

-- one row for each unique list of tables
CREATE TABLE val.lists
  (lid SERIAL, 
  name TEXT NOT NULL, 
  comment TEXT NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT lists_pk PRIMARY KEY (lid), 
  CONSTRAINT lists_unique_name   UNIQUE (name) ); 
GRANT SELECT ON val.lists TO PUBLIC;
GRANT INSERT ON val.lists TO manager_role;
GRANT UPDATE ON val.lists_lid_seq TO manager_role;

-- relationship between a list of tables (val.lists) and tables (val.tables)
CREATE TABLE val.tablelists
  (lid INT, 
  tid INT, 
  CONSTRAINT tablelists_pk PRIMARY KEY (lid,tid),
  CONSTRAINT tablelists_lid_fk FOREIGN KEY (lid) REFERENCES val.lists(lid),
  CONSTRAINT tablelists_tid_fk FOREIGN KEY (tid) REFERENCES val.tables(tid) );
GRANT SELECT ON val.tablelists TO PUBLIC;
GRANT INSERT ON val.tablelists TO manager_role;

-- each purpose may have several versions
CREATE TABLE val.versions
  (vid SERIAL, 
  pid INTEGER NOT NULL, 
  lid INTEGER NOT NULL, 
  major INTEGER NOT NULL, 
  minor INTEGER NOT NULL, 
  comment TEXT NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT vid_pk PRIMARY KEY (vid),
  CONSTRAINT vid_pid_fk FOREIGN KEY (pid) REFERENCES val.purposes(pid),
  CONSTRAINT vid_lid_fk FOREIGN KEY (lid) REFERENCES val.lists(lid), 
  CONSTRAINT versions_unique_combo   UNIQUE (pid,major,minor) );
GRANT SELECT ON val.versions TO PUBLIC;
GRANT INSERT ON val.versions TO manager_role;
GRANT UPDATE ON val.versions_vid_seq TO manager_role;

-- each purpose/extension can be extended
CREATE TABLE val.extensions
  (eid SERIAL, 
  vid INTEGER NOT NULL, 
  extension INTEGER NOT NULL,
  create_time TIMESTAMP WITH TIME ZONE NOT NULL, 
  create_user TEXT NOT NULL,  
  CONSTRAINT extensions_pk PRIMARY KEY (eid),
  CONSTRAINT extensions_vid_fk FOREIGN KEY (vid) REFERENCES val.versions(vid),
  CONSTRAINT extensions_unique_two UNIQUE (vid,extension) );
GRANT SELECT ON val.extensions TO PUBLIC;
GRANT INSERT ON val.extensions TO manager_role;
GRANT UPDATE ON val.extensions_eid_seq TO manager_role;

-- for each extension, a list of groups in the extension
CREATE TABLE val.extensionlists
  (eid INT, 
  gid INT, 
  CONSTRAINT extensionlists_pk PRIMARY KEY (eid,gid),
  CONSTRAINT extensionlists_lid_fk FOREIGN KEY (eid) REFERENCES val.extensions(eid),
  CONSTRAINT extensionlists_tid_fk FOREIGN KEY (gid) REFERENCES val.groups(gid) );
GRANT SELECT ON val.extensionlists TO PUBLIC;
GRANT INSERT ON val.extensionlists TO manager_role;

--
-- tst schema tables
--

CREATE TABLE tst.calib1 
  (cid INTEGER, channel INTEGER, flag INTEGER , dtoe NUMERIC,
   CONSTRAINT tst_calib1_pk PRIMARY KEY (cid,channel) );
GRANT SELECT ON tst.calib1 TO PUBLIC;
GRANT INSERT ON tst.calib1 TO val_role;

CREATE TABLE tst.calib2 
  (cid INTEGER, channel INTEGER, status TEXT, 
   CONSTRAINT tst_calib2_pk PRIMARY KEY (cid,channel) );
GRANT SELECT ON tst.calib2 TO PUBLIC;
GRANT INSERT ON tst.calib2 TO val_role;

--
-- trk schema tables
--

-- initial calibration tables

CREATE TABLE trk.preampstraw
  (cid INTEGER, index INTEGER,  
   delay_hv NUMERIC, delay_cal NUMERIC,  
   threshold_hv NUMERIC, threshold_cal NUMERIC, gain NUMERIC,
   CONSTRAINT trk_preampstraw_pk PRIMARY KEY (cid,index) );
GRANT SELECT ON trk.preampstraw TO PUBLIC;
GRANT INSERT ON trk.preampstraw TO trk_role;


CREATE TABLE trk.preamprstraw
  (cid INTEGER, index INTEGER,  
   delay_hv NUMERIC, delay_cal NUMERIC,  
   threshold_hv NUMERIC, threshold_cal NUMERIC, gain NUMERIC,
   CONSTRAINT trk_preamprstraw_pk PRIMARY KEY (cid,index) );
GRANT SELECT ON trk.preamprstraw TO PUBLIC;
GRANT INSERT ON trk.preamprstraw TO trk_role;

CREATE TABLE trk.delaypanel
  (cid INTEGER, index INTEGER, delay NUMERIC,
   CONSTRAINT trk_delaypanel_pk PRIMARY KEY (cid,index) );
GRANT SELECT ON trk.delaypanel TO PUBLIC;
GRANT INSERT ON trk.delaypanel TO trk_role;


CREATE TABLE trk.thresholdrstraw
  (cid INTEGER, index INTEGER, 
   threshold_hv NUMERIC, threshold_cal NUMERIC,
   CONSTRAINT trk_thresholdrstraw_pk PRIMARY KEY (cid,index) );
GRANT SELECT ON trk.thresholdrstraw TO PUBLIC;
GRANT INSERT ON trk.thresholdrstraw TO trk_role;

-- tracker alignment

CREATE TABLE trk.aligntracker
  (cid INTEGER, index INTEGER, strawid TEXT, 
   dx NUMERIC, dy NUMERIC, dz NUMERIC, 
   rx NUMERIC, ry NUMERIC, rz NUMERIC, 
   CONSTRAINT trk_aligntracker_pk PRIMARY KEY (cid,index) );
GRANT SELECT ON trk.aligntracker TO PUBLIC;
GRANT INSERT ON trk.aligntracker TO trk_role;

CREATE TABLE trk.alignplane
  (cid INTEGER, index INTEGER, strawid TEXT, 
   dx NUMERIC, dy NUMERIC, dz NUMERIC, 
   rx NUMERIC, ry NUMERIC, rz NUMERIC, 
   CONSTRAINT trk_alignplane_pk PRIMARY KEY (cid,index) );
GRANT SELECT ON trk.alignplane TO PUBLIC;
GRANT INSERT ON trk.alignplane TO trk_role;

CREATE TABLE trk.alignpanel
  (cid INTEGER, index INTEGER, strawid TEXT, 
   dx NUMERIC, dy NUMERIC, dz NUMERIC, 
   rx NUMERIC, ry NUMERIC, rz NUMERIC, 
   CONSTRAINT trk_alignpanel_pk PRIMARY KEY (cid,index) );
GRANT SELECT ON trk.alignpanel TO PUBLIC;
GRANT INSERT ON trk.alignpanel TO trk_role;

CREATE TABLE trk.alignstraw
  (cid INTEGER, index INTEGER, StrawId TEXT, 
   wire_cal_dV NUMERIC, wire_cal_dW NUMERIC,
   wire_hv_dV NUMERIC, wire_hv_dW NUMERIC,
   straw_cal_dV NUMERIC, straw_cal_dW NUMERIC,
   straw_hv_dV NUMERIC, straw_hv_dW NUMERIC,
   CONSTRAINT trk_alignstraw_pk PRIMARY KEY (cid,index) );
GRANT SELECT ON trk.alignstraw TO PUBLIC;
GRANT INSERT ON trk.alignstraw TO trk_role;

-- tracker component status

CREATE TABLE trk.panelstatus
  (cid INTEGER, strawid TEXT, strawstatus TEXT,
   CONSTRAINT trk_panelstatus_pk PRIMARY KEY (cid,strawid) );
GRANT SELECT ON trk.panelstatus TO PUBLIC;
GRANT INSERT ON trk.panelstatus TO trk_role;

CREATE TABLE trk.planestatus
  (cid INTEGER, strawid TEXT, strawstatus TEXT,
   CONSTRAINT trk_planestatus_pk PRIMARY KEY (cid,strawid) );
GRANT SELECT ON trk.planestatus TO PUBLIC;
GRANT INSERT ON trk.planestatus TO trk_role;

CREATE TABLE trk.strawstatusshort
  (cid INTEGER, strawid TEXT, strawstatus TEXT,
   CONSTRAINT trk_strawstatusshort_pk PRIMARY KEY (cid,strawid) );
GRANT SELECT ON trk.strawstatusshort TO PUBLIC;
GRANT INSERT ON trk.strawstatusshort TO trk_role;

CREATE TABLE trk.strawstatuslong
  (cid INTEGER, strawid TEXT, strawstatus TEXT,
   CONSTRAINT trk_strawstatuslong_pk PRIMARY KEY (cid,strawid) );
GRANT SELECT ON trk.strawstatuslong TO PUBLIC;
GRANT INSERT ON trk.strawstatuslong TO trk_role;

-- tracker calibration table added 7/2021

CREATE TABLE trk.delayrstraw
  (cid INTEGER, 
   straw INTEGER,  delay_hv NUMERIC, delay_cal NUMERIC,  
   CONSTRAINT trk_delayrstraw_pk PRIMARY KEY (cid,straw) );
GRANT SELECT ON trk.delayrstraw TO PUBLIC;
GRANT INSERT ON trk.delayrstraw TO trk_role;

--
-- sim schema tables
--

-- table of mixing efficiencies

CREATE TABLE sim.efficiencies
  (cid INTEGER, 
   tag TEXT, numerator INTEGER, denominator INTEGER, eff NUMERIC,  
   CONSTRAINT sim_efficiencies_pk PRIMARY KEY (cid,tag) );
GRANT SELECT ON sim.efficiencies TO PUBLIC;
GRANT INSERT ON sim.efficiencies TO sim_role;

CREATE TABLE sim.efficiencies2
  (cid INTEGER, 
   tag TEXT, numerator BIGINT, denominator BIGINT, eff NUMERIC,  
   CONSTRAINT sim_efficiencies2_pk PRIMARY KEY (cid,tag) );
GRANT SELECT ON sim.efficiencies2 TO PUBLIC;
GRANT INSERT ON sim.efficiencies2 TO sim_role;
