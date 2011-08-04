#pragma once
namespace swarm {

	// Note: the header _MUST_ be padded to 16-byte boundary
	struct ALIGN(16) swarm_header
	{
		char magic[6];		// Magic string to quickly verify this is a swarm file (== 'SWARM\0')
		char version[2];	// File format version
		char m_type[76];	// user-defined file content ID/description
		uint32_t flags;		// user-defined flags
		uint64_t datalen;	// length of data in the file (0xFFFFFFFFFFFFFFFF for unknown)

		static const uint64_t npos = 0xFFFFFFFFFFFFFFFFLL;

		swarm_header(const std::string &type, int flags_ = 0, uint64_t datalen_ = npos)
		{
			strcpy(magic, "SWARM");
			strcpy(version, "0");

			strncpy(this->m_type, type.c_str(), sizeof(this->m_type));
			this->m_type[sizeof(this->m_type)-1] = '\0';

			flags = flags_;
			datalen = datalen_;
		}

		std::string type() const
		{
			const char *c = strstr(m_type, "//");
			if(!c) { return trim(m_type); }
			return trim(std::string(m_type, c - m_type));
		}

		bool is_compatible(const swarm_header &a)
		{
			bool ffver_ok = memcmp(magic, a.magic, 8) == 0;
			std::string t = type();
			bool type_ok = t.empty() || (t == a.type());
			return ffver_ok && type_ok;
		}
	};

	struct ALIGN(16) swarm_index_header : public swarm_header
	{
		uint64_t	timestamp;	// datafile timestamp (mtime)
		uint64_t	datafile_size;	// datafile file size

		swarm_index_header(const std::string &type, uint64_t timestamp_ = 0, uint64_t datafile_size_ = 0 , uint64_t datalen_ = npos)
			: swarm_header(type, 0, datalen_), timestamp(timestamp_), datafile_size(datafile_size_)
		{
		}
	};

}
