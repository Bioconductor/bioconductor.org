#!/usr/bin/env ruby

require 'rubygems'
require 'right_aws'
require 'pp'
require 'yaml'


path = File.expand_path $0
segs = path.split "/"
segs.pop
path = segs.join "/"

credentials = YAML.load_file("#{path}/credentials.yaml")
access_key = credentials["access_key"]
secret_key = credentials["secret_key"]


ec2   = RightAws::Ec2.new(access_key, secret_key)
amis = ec2.describe_images_by_owner('self')
me = amis.first[:aws_owner]


never_delete_these_snapshots = ["snap-21e2964b"]

ami_ids = amis.map{|i|i[:aws_id]}

#instances = ec2.describe_instances

snaps_to_keep = never_delete_these_snapshots

for ami in amis
  for bdm in ami[:block_device_mappings]
    snaps_to_keep.push bdm[:ebs_snapshot_id]
  end
end

snaps_to_keep.sort.uniq!

all_snaps = ec2.describe_snapshots

snaps_to_delete = []
for snap in all_snaps
  unless snaps_to_keep.include? snap[:aws_id]
    if snap[:aws_owner] == me # ??????????
      snaps_to_delete.push(snap)
    end
  end
end

snap_ids_to_delete = snaps_to_delete.map{|i|i[:aws_id]}


for id in snap_ids_to_delete
  begin
    ec2.delete_snapshot id
    puts "deleted snapshot #{id}"
  rescue RightAws::AwsError => ex
    puts "error deleting snapshot #{id}"
  end
end

volumes = ec2.describe_volumes


volumes_to_delete = []

for vol in volumes
  next if vol[:aws_attachment_status] == "attached"
  next if vol[:aws_status] == "in-use"
  # think about whether this should run or not:
  ######next if snaps_to_keep.include? vol[:snapshot_id]
  volumes_to_delete.push vol
end

volume_ids_to_delete = volumes_to_delete.map{|i|i[:aws_id]}

puts volume_ids_to_delete.size


for vol in volume_ids_to_delete
  res = ec2.delete_volume(vol)
  puts "result of deleting #{vol}: #{res}"
end
