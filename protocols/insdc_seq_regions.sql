-- Provided by James Allen Nov2018

-- Labelling contigs as 'ENA' enables link outs in the web display
insert into seq_region_attrib
  select seq_region_id, 317, "ENA" from
    seq_region inner join
    coord_system using (coord_system_id)
  where coord_system.name = "contig";

-- If we don't use the INSDC name as the seq_region name, adding it
-- as a synonym makes it searchable. This approach assumes that the
-- seq_region name is nonetheless present in the INSDC record.
-- Need to download the mapping file(s) from NCBI and save it locally.
create temporary table tmp_map (acc varchar(255), id varchar(255));
load data local infile "scaffold_localID2acc" into table tmp_map;
insert into seq_region_synonym (seq_region_id, synonym, external_db_id)
  select seq_region_id, id, 50710 from
    seq_region_synonym inner join
    tmp_map on synonym = acc;
drop temporary table tmp_map;
