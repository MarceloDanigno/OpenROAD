#BSD 3-Clause License
#
#Copyright (c) 2023, The Regents of the University of California
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
#3. Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

sta::define_cmd_args "sec_evaluate" {}

proc sec_evaluate { args } {

  sec::evaluate

}

sta::define_cmd_args "sec_read_design_metrics" {
    [-json filename]
}

proc sec_read_design_metrics { args } {
  sta::parse_key_args "sec_read_design_metrics" args keys \
      { -json }

  set json_file ""
  if { [info exists keys(-json)] } {
    set json_file $keys(-json)
  }

  sec::read_design_metrics $json_file

}


sta::define_cmd_args "sec_read_baseline_metrics" {
    [-json filename]
}

proc sec_read_baseline_metrics { args } {
  sta::parse_key_args "sec_read_baseline_metrics" args keys \
      { -json }

  set json_file ""
  if { [info exists keys(-json)] } {
    set json_file $keys(-json)
  }

  sec::read_baseline_metrics $json_file

}

sta::define_cmd_args "sec_compute_overall" {}

proc sec_compute_overall { args } {

  sec::compute_overall

}
